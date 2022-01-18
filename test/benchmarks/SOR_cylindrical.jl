using BenchmarkTools
using SolidStateDetectors
# using Plots
using Unitful

using SolidStateDetectors: rb_even, rb_odd

### Inner loop of SOR 
T = Float32
# sim = Simulation{T}(SSD_examples[:InvertedCoax])
sim = Simulation{T}(SSD_examples[:Coax])

# imp_density = SolidStateDetectors.CylindricalImpurityDensity{Float32}((0.0f0, 0.0f0, 1f16), (0.0f0, 0.0f0, 0.0f0));
# sim.detector = SolidStateDetector(sim.detector, imp_density); 
# sim.detector = SolidStateDetector(sim.detector, contact_id = 1, contact_potential = 2000); # Make it undepleted

calculate_electric_potential!(sim, depletion_handling = false,
    convergence_limit = 1e-3,
    use_nthreads = 4,
    max_tick_distance = (2u"mm", 1u"°", 2u"mm"),
    # refinement_limits = missing,
    refinement_limits = [0.2, 0.1, 0.05, 0.02, 0.01],
)
# plot(
#     plot(sim.point_types, φ = 30),
#     plot(sim.electric_potential, φ = 30),
# )

pssrb = SolidStateDetectors.PotentialSimulationSetupRB(sim.detector,
    sim.electric_potential.grid,
    SolidStateDetectors.material_properties[SolidStateDetectors.materials["vacuum"]],
    sim.electric_potential.data);


nthreads = 1 #Base.Threads.nthreads()
is_weighting_potential = Val{false}()
depletion_handling = Val{false}()
only2d = Val{false}()
update_even_points = Val{true}()
update_uneven_points = Val{false}()

begin
    @info "Grid size: $(size(sim.electric_potential.data))"
    # SolidStateDetectors.outerloop!(pssrb, nthreads, update_even_points, depletion_handling, is_weighting_potential, only2d)
    for dp in (false, true)
        depletion_handling = Val{dp}()
        @info "Depletion handling: $dp"
        for nt in (1, 2, 4, 8, 16, 32, 64)
            @info "N Threads: $nt"
            @btime SolidStateDetectors.outerloop!($pssrb, $nt, $update_even_points, $depletion_handling, $is_weighting_potential, $only2d)
        end
    end
end


# middleloop!
idx3 = 2
even_points = true
update_even_points = Val(even_points)
is_weighting_potential = Val(false)
rb_tar_idx, rb_src_idx = even_points ? (rb_even, rb_odd) : (rb_odd, rb_even)
depletion_handling = Val{false}()
idx3iseven = Val(iseven(idx3))

SolidStateDetectors.middleloop!(idx3, rb_tar_idx, rb_src_idx, pssrb,
    update_even_points, depletion_handling, is_weighting_potential, only2d, idx3iseven)

@code_warntype SolidStateDetectors.middleloop!(idx3, rb_tar_idx, rb_src_idx, pssrb,
    update_even_points, depletion_handling, is_weighting_potential, only2d, idx3iseven)

@btime SolidStateDetectors.middleloop!($idx3, $rb_tar_idx, $rb_src_idx, $pssrb,
    $update_even_points, $depletion_handling, $is_weighting_potential, $only2d, $idx3iseven)



# Innerloop:
line_weights = Array{T,2}(undef, size(pssrb.potential, 1) - 2, 6)

depletion_handling = Val{true}()
idx3_is_even = Val(iseven(idx3))

ir = idx3
inr = ir - 1

pwwrr = pssrb.geom_weights[1][1, inr]
pwwrl = pssrb.geom_weights[1][2, inr]
r_inv_pwΔmpr = pssrb.geom_weights[1][3, inr]
Δr_ext_inv_r_pwmprr = pssrb.geom_weights[1][4, inr]
Δr_ext_inv_l_pwmprl = pssrb.geom_weights[1][5, inr]
Δmpr_squared = pssrb.geom_weights[1][6, inr]

iφ = 2
inφ = iφ - 1

pwwφr = pssrb.geom_weights[2][1, inφ]
pwwφl = pssrb.geom_weights[2][2, inφ]
pwΔmpφ = pssrb.geom_weights[2][3, inφ]
Δφ_ext_inv_r = pssrb.geom_weights[2][4, iφ]
Δφ_ext_inv_l = pssrb.geom_weights[2][4, inφ]

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

SolidStateDetectors.load_weights_for_innerloop!(
    line_weights,
    pssrb, iφ, inφ, ir, inr,
    update_even_points, idx3_is_even,
    pwwrr, pwwrl, pwwφr, pwwφl,
    pwwrr_pwwφr, pwwrl_pwwφr, pwwrr_pwwφl, pwwrl_pwwφl,
    pwΔmpφ_Δmpr_squared,
    Δr_ext_inv_r_pwmprr_pwΔmpφ, Δr_ext_inv_l_pwmprl_pwΔmpφ,
    r_inv_pwΔmpr_Δφ_ext_inv_r, r_inv_pwΔmpr_Δφ_ext_inv_l,
)

@code_warntype SolidStateDetectors.load_weights_for_innerloop!(
    line_weights,
    pssrb, iφ, inφ, ir, inr,
    update_even_points, idx3_is_even,
    pwwrr, pwwrl, pwwφr, pwwφl,
    pwwrr_pwwφr, pwwrl_pwwφr, pwwrr_pwwφl, pwwrl_pwwφl,
    pwΔmpφ_Δmpr_squared,
    Δr_ext_inv_r_pwmprr_pwΔmpφ, Δr_ext_inv_l_pwmprl_pwΔmpφ,
    r_inv_pwΔmpr_Δφ_ext_inv_r, r_inv_pwΔmpr_Δφ_ext_inv_l,
)

@btime SolidStateDetectors.load_weights_for_innerloop!(
    $line_weights,
    $pssrb, $iφ, $inφ, $ir, $inr,
    $update_even_points, $idx3_is_even,
    $pwwrr, $pwwrl, $pwwφr, $pwwφl,
    $pwwrr_pwwφr, $pwwrl_pwwφr, $pwwrr_pwwφl, $pwwrl_pwwφl,
    $pwΔmpφ_Δmpr_squared,
    $Δr_ext_inv_r_pwmprr_pwΔmpφ, $Δr_ext_inv_l_pwmprl_pwΔmpφ,
    $r_inv_pwΔmpr_Δφ_ext_inv_r, $r_inv_pwΔmpr_Δφ_ext_inv_l,
)

SolidStateDetectors.innerloop!(line_weights, pssrb, iφ, inφ, ir, inr, rb_tar_idx, rb_src_idx,
    update_even_points, idx3_is_even,
    depletion_handling, is_weighting_potential, only2d)

@code_warntype SolidStateDetectors.innerloop!(line_weights, pssrb, iφ, inφ, ir, inr, rb_tar_idx, rb_src_idx,
    update_even_points, idx3_is_even,
    depletion_handling, is_weighting_potential, only2d)

@btime SolidStateDetectors.innerloop!($line_weights, $pssrb, $iφ, $inφ, $ir, $inr, $rb_tar_idx, $rb_src_idx,
    $update_even_points, $idx3_is_even,
    $depletion_handling, $is_weighting_potential, $only2d)

@code_llvm SolidStateDetectors.innerloop!(line_weights, pssrb, iφ, inφ, ir, inr, rb_tar_idx, rb_src_idx,
    update_even_points, idx3_is_even,
    depletion_handling, is_weighting_potential, only2d)



# apply_boundary_conditions!
begin
    @btime SolidStateDetectors.apply_boundary_conditions!($pssrb, $update_even_points, $only2d)
    @btime SolidStateDetectors.apply_boundary_conditions!($pssrb, $update_uneven_points, $only2d)
end
