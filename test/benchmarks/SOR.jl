using BenchmarkTools
using SolidStateDetectors
using Plots
using Unitful

using SolidStateDetectors: rb_even, rb_odd

### Inner loop of SOR 
begin
    T = Float32
    sim = Simulation{T}(SSD_examples[:Coax]);

    imp_density = SolidStateDetectors.CylindricalImpurityDensity{Float32}((0.0f0, 0.0f0, 1f16), (0.0f0, 0.0f0, 0.0f0));
    sim.detector = SolidStateDetector(sim.detector, imp_density); 
    sim.detector = SolidStateDetector(sim.detector, contact_id = 1, contact_potential = 2000); # Make it undepleted

    calculate_electric_potential!( sim, depletion_handling = false,
                                        convergence_limit = 1e-4,
                                        max_tick_distance = 1u"mm",
                                        # refinement_limits = [0.2, 0.1, 0.05],
                                        )

    plot(
        plot(sim.point_types, φ = 30),
        plot(sim.electric_potential, φ = 30),
    )

    pssrb = SolidStateDetectors.PotentialSimulationSetupRB( sim.detector, 
                                                            sim.electric_potential.grid, 
                                                            # Grid(sim, max_tick_distance = 1u"mm"), 
                                                            SolidStateDetectors.material_properties[SolidStateDetectors.materials["vacuum"]],
                                                            sim.electric_potential.data, 
                                                            sor_consts = (T(1), T(1)) );


    nthreads = 1 #Base.Threads.nthreads()
    is_weighting_potential = Val{false}()
    depletion_handling = Val{false}()
    only2d = Val{false}()
    update_even_points = Val{true}()
    update_uneven_points = Val{false}()
end

begin
    @info size(sim.electric_potential.data)
    SolidStateDetectors.outerloop!(pssrb, nthreads, update_even_points, depletion_handling, is_weighting_potential, only2d)
    depletion_handling = Val{false}()
    # @code_warntype SolidStateDetectors.outerloop!(pssrb, nthreads, update_even_points, depletion_handling, is_weighting_potential, only2d)
    nthreads = 1 #Base.Threads.nthreads()
    @btime SolidStateDetectors.outerloop!($pssrb, $nthreads, $update_even_points, $depletion_handling, $is_weighting_potential, $only2d)
    nthreads = 4 #Base.Threads.nthreads()
    @btime SolidStateDetectors.outerloop!($pssrb, $nthreads, $update_even_points, $depletion_handling, $is_weighting_potential, $only2d)
    # 164 μs v0.6.0
end
begin
    depletion_handling = Val{true}()
    @btime SolidStateDetectors.outerloop!($pssrb, $nthreads, $update_even_points, $depletion_handling, $is_weighting_potential, $only2d)
end

# middleloop!
begin 
    idx3 = 2
    even_points = true
    update_even_points = Val(even_points)
    is_weighting_potential = Val(false)
    rb_tar_idx, rb_src_idx = even_points ? (rb_even, rb_odd) : (rb_odd, rb_even)
    gw1 = pssrb.geom_weights[1].weights;  # r or x 
    gw2 = pssrb.geom_weights[2].weights;  # φ or y
    gw3 = pssrb.geom_weights[3].weights;  # z or z
    depletion_handling = Val{false}()
    idx3iseven = Val(iseven(idx3))

    SolidStateDetectors.middleloop!(idx3, rb_tar_idx, rb_src_idx, gw1, gw2, gw3, pssrb, 
                         update_even_points, depletion_handling, is_weighting_potential, only2d, idx3iseven)
    
    @code_warntype SolidStateDetectors.middleloop!(idx3, rb_tar_idx, rb_src_idx, gw1, gw2, gw3, pssrb, 
                         update_even_points, depletion_handling, is_weighting_potential, only2d, idx3iseven)
    

    @btime SolidStateDetectors.middleloop!($idx3, $rb_tar_idx, $rb_src_idx, $gw1, $gw2, $gw3, $pssrb, 
                         $update_even_points, $depletion_handling, $is_weighting_potential, $only2d, $idx3iseven)

end



# Innerloop:
begin
    depletion_handling = Val{false}()
    idx3_is_even = Val(iseven(idx3))
    ir = idx3
    inr = ir - 1
    gw_r = gw1;
    pwwrr               = gw_r[1, inr]
    pwwrl               = gw_r[2, inr]
    r_inv_pwΔmpr        = gw_r[3, inr]
    Δr_ext_inv_r_pwmprr = gw_r[4, inr] 
    Δr_ext_inv_l_pwmprl = gw_r[5, inr] 
    Δmpr_squared        = gw_r[6, inr]  
    iφ = 2
    inφ = iφ - 1
    gw_φ = gw2;
    pwwφr      = gw_φ[1, inφ]
    pwwφl        = gw_φ[2, inφ]
    pwΔmpφ       = gw_φ[3, inφ]
    Δφ_ext_inv_r = gw_φ[4,  iφ]
    Δφ_ext_inv_l = gw_φ[4, inφ]
    if inr == 1
        pwwφr = T(0.5)
        pwwφl = T(0.5)
        pwΔmpφ = T(2π)
        Δφ_ext_inv_r = inv(pwΔmpφ)
        Δφ_ext_inv_l = Δφ_ext_inv_r
    end
    pwwrr_pwwφr = pwwrr * pwwφr
    pwwrr_pwwφl = pwwrr * pwwφl
    pwwrl_pwwφr = pwwrl * pwwφr
    pwwrl_pwwφl = pwwrl * pwwφl

    Δr_ext_inv_r_pwmprr_pwΔmpφ = Δr_ext_inv_r_pwmprr * pwΔmpφ
    Δr_ext_inv_l_pwmprl_pwΔmpφ = Δr_ext_inv_l_pwmprl * pwΔmpφ
    pwΔmpφ_Δmpr_squared = pwΔmpφ * Δmpr_squared
    r_inv_pwΔmpr_Δφ_ext_inv_r = r_inv_pwΔmpr * Δφ_ext_inv_r
    r_inv_pwΔmpr_Δφ_ext_inv_l = r_inv_pwΔmpr * Δφ_ext_inv_l
    gw_z = gw3;

    SolidStateDetectors.load_weights!(pssrb, rb_tar_idx, rb_src_idx, gw_z, ir, inr, iφ, inφ, 
                pwwφr, pwwφl, pwwrl, pwwrr,
                update_even_points, even_points, idx3_is_even, 
                Δr_ext_inv_r_pwmprr_pwΔmpφ, Δr_ext_inv_l_pwmprl_pwΔmpφ, pwΔmpφ_Δmpr_squared, 
                r_inv_pwΔmpr_Δφ_ext_inv_r, r_inv_pwΔmpr_Δφ_ext_inv_l,
                pwwrr_pwwφr, pwwrr_pwwφl, pwwrl_pwwφr, pwwrl_pwwφl, 
                depletion_handling, is_weighting_potential, only2d)
    @btime SolidStateDetectors.load_weights!($pssrb, $rb_tar_idx, $rb_src_idx, $gw_z, $ir, $inr, $iφ, $inφ, 
                $pwwφr, $pwwφl, $pwwrl, $pwwrr,
                $update_even_points, $even_points, $idx3_is_even, 
                $Δr_ext_inv_r_pwmprr_pwΔmpφ, $Δr_ext_inv_l_pwmprl_pwΔmpφ, $pwΔmpφ_Δmpr_squared, 
                $r_inv_pwΔmpr_Δφ_ext_inv_r, $r_inv_pwΔmpr_Δφ_ext_inv_l,
                $pwwrr_pwwφr, $pwwrr_pwwφl, $pwwrl_pwwφr, $pwwrl_pwwφl, 
                $depletion_handling, $is_weighting_potential, $only2d)

    SolidStateDetectors.innerloop!(pssrb, rb_tar_idx, rb_src_idx, gw_z, ir, inr, iφ, inφ, 
                pwwφr, pwwφl, pwwrl, pwwrr,
                update_even_points, even_points, idx3_is_even, 
                Δr_ext_inv_r_pwmprr_pwΔmpφ, Δr_ext_inv_l_pwmprl_pwΔmpφ, pwΔmpφ_Δmpr_squared, 
                r_inv_pwΔmpr_Δφ_ext_inv_r, r_inv_pwΔmpr_Δφ_ext_inv_l,
                pwwrr_pwwφr, pwwrr_pwwφl, pwwrl_pwwφr, pwwrl_pwwφl, 
                depletion_handling, is_weighting_potential, only2d)
    
    @code_warntype SolidStateDetectors.innerloop!(pssrb, rb_tar_idx, rb_src_idx, gw_z, ir, inr, iφ, inφ,
                pwwφr, pwwφl, pwwrl, pwwrr,
                update_even_points, even_points, idx3_is_even, 
                Δr_ext_inv_r_pwmprr_pwΔmpφ, Δr_ext_inv_l_pwmprl_pwΔmpφ, pwΔmpφ_Δmpr_squared, 
                r_inv_pwΔmpr_Δφ_ext_inv_r, r_inv_pwΔmpr_Δφ_ext_inv_l,
                pwwrr_pwwφr, pwwrr_pwwφl, pwwrl_pwwφr, pwwrl_pwwφl, 
                depletion_handling, is_weighting_potential, only2d)

    @btime SolidStateDetectors.innerloop!($pssrb, $rb_tar_idx, $rb_src_idx, $gw_z, $ir, $inr, $iφ, $inφ, 
                $pwwφr, $pwwφl, $pwwrl, $pwwrr,
                $update_even_points, $even_points, $idx3_is_even, 
                $Δr_ext_inv_r_pwmprr_pwΔmpφ, $Δr_ext_inv_l_pwmprl_pwΔmpφ, $pwΔmpφ_Δmpr_squared, 
                $r_inv_pwΔmpr_Δφ_ext_inv_r, $r_inv_pwΔmpr_Δφ_ext_inv_l,
                $pwwrr_pwwφr, $pwwrr_pwwφl, $pwwrl_pwwφr, $pwwrl_pwwφl, 
                $depletion_handling, $is_weighting_potential, $only2d)


    @code_llvm SolidStateDetectors.innerloop!(pssrb, rb_tar_idx, rb_src_idx, gw_z, ir, inr, iφ, inφ,
                pwwφr, pwwφl, pwwrl, pwwrr,
                update_even_points, even_points, idx3_is_even, 
                Δr_ext_inv_r_pwmprr_pwΔmpφ, Δr_ext_inv_l_pwmprl_pwΔmpφ, pwΔmpφ_Δmpr_squared, 
                r_inv_pwΔmpr_Δφ_ext_inv_r, r_inv_pwΔmpr_Δφ_ext_inv_l,
                pwwrr_pwwφr, pwwrr_pwwφl, pwwrl_pwwφr, pwwrl_pwwφl, 
                depletion_handling, is_weighting_potential, only2d)
end

# apply_boundary_conditions!
begin 
    @btime SolidStateDetectors.apply_boundary_conditions!($pssrb, $update_even_points, $only2d)
    @btime SolidStateDetectors.apply_boundary_conditions!($pssrb, $update_uneven_points, $only2d)
end

begin
    depletion_handling = Val{false}()
    SolidStateDetectors.update!(pssrb; use_nthreads = nthreads, depletion_handling = depletion_handling,
        only2d = only2d, is_weighting_potential= is_weighting_potential)
    @btime SolidStateDetectors.update!($pssrb; use_nthreads = $nthreads, depletion_handling = $depletion_handling,
        only2d = $only2d, is_weighting_potential= $is_weighting_potential)
  
    depletion_handling = Val{true}()
    SolidStateDetectors.update!(pssrb; use_nthreads = nthreads, depletion_handling = depletion_handling,
        only2d = only2d, is_weighting_potential= is_weighting_potential)
    @btime SolidStateDetectors.update!($pssrb; use_nthreads = $nthreads, depletion_handling = $depletion_handling,
        only2d = $only2d, is_weighting_potential= $is_weighting_potential)
end