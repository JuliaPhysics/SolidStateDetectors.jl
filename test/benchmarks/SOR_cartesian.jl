using BenchmarkTools
using SolidStateDetectors
using Plots
using Unitful

using SolidStateDetectors: rb_even, rb_odd

### Inner loop of SOR 
T = Float32
sim = Simulation{T}(SSD_examples[:CGD]);

calculate_electric_potential!( sim, depletion_handling = false,
                                    convergence_limit = 1e-4,
                                    max_tick_distance = 1u"mm",
                                    missing,
                                    # refinement_limits = [0.2, 0.1, 0.05, 0.02],
                                    refinement_limits = missing, 
                                    )

# plot(
#     plot(sim.point_types, y = 0),
#     plot(sim.electric_potential, y = 0),
# )

pssrb = SolidStateDetectors.PotentialSimulationSetupRB( sim.detector, 
                                                        sim.electric_potential.grid, 
                                                        # Grid(sim, max_tick_distance = 1u"mm"), 
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



# load_weights_for_innerloop! & innerloop! 
depletion_handling = Val{false}()

idx3_is_even = Val(iseven(idx3))

iz = idx3
inz = iz - 1 
            
pwwzr        = pssrb.geom_weights[3][1, inz]
pwwzl        = pssrb.geom_weights[3][2, inz]
pwΔmpz       = pssrb.geom_weights[3][3, inz]
Δz_ext_inv_r = pssrb.geom_weights[3][4, inz + 1]
Δz_ext_inv_l = pssrb.geom_weights[3][4, inz]

line_weights = Array{T, 2}(undef, size(pssrb.potential, 1) - 2, 6)

iy = 2
iny = iy - 1

pwwyr  = pssrb.geom_weights[2][1, iny]
pwwyl  = pssrb.geom_weights[2][2, iny]
pwΔmpy = pssrb.geom_weights[2][3, iny] 
pwΔmpy_pwΔmpz = pwΔmpy * pwΔmpz
Δy_ext_inv_r_pwΔmpz  = pssrb.geom_weights[2][4, iny + 1] * pwΔmpz
Δy_ext_inv_l_pwΔmpz  = pssrb.geom_weights[2][4, iny]     * pwΔmpz
Δz_ext_inv_r_pwΔmpy = Δz_ext_inv_r * pwΔmpy
Δz_ext_inv_l_pwΔmpy = Δz_ext_inv_l * pwΔmpy

pwwyr_pwwzr = pwwyr * pwwzr
pwwyr_pwwzl = pwwyr * pwwzl
pwwyl_pwwzr = pwwyl * pwwzr
pwwyl_pwwzl = pwwyl * pwwzl

SolidStateDetectors.load_weights_for_innerloop!(line_weights, pssrb, iy, iny, iz, inz,
            update_even_points, idx3_is_even, 
            pwwzr, pwwzl, pwwyr, pwwyl,
            pwwyr_pwwzr, pwwyr_pwwzl, pwwyl_pwwzr, pwwyl_pwwzl,
            pwΔmpy_pwΔmpz,
            Δy_ext_inv_r_pwΔmpz, Δy_ext_inv_l_pwΔmpz, Δz_ext_inv_r_pwΔmpy, Δz_ext_inv_l_pwΔmpy)


@btime SolidStateDetectors.load_weights_for_innerloop!($line_weights, $pssrb, $iy, $iny, $iz, $inz,
            $update_even_points, $idx3_is_even, 
            $pwwzr, $pwwzl, $pwwyr, $pwwyl,
            $pwwyr_pwwzr, $pwwyr_pwwzl, $pwwyl_pwwzr, $pwwyl_pwwzl,
            $pwΔmpy_pwΔmpz,
            $Δy_ext_inv_r_pwΔmpz, $Δy_ext_inv_l_pwΔmpz, $Δz_ext_inv_r_pwΔmpy, $Δz_ext_inv_l_pwΔmpy)

SolidStateDetectors.innerloop!(line_weights, pssrb, iy, iny, iz, inz, rb_tar_idx, rb_src_idx,
            update_even_points, idx3_is_even, 
            depletion_handling, is_weighting_potential, only2d)

@code_warntype SolidStateDetectors.innerloop!(line_weights, pssrb, iy, iny, iz, inz, rb_tar_idx, rb_src_idx,
            update_even_points, idx3_is_even, 
            depletion_handling, is_weighting_potential, only2d)

@btime SolidStateDetectors.innerloop!($line_weights, $pssrb, $iy, $iny, $iz, $inz, $rb_tar_idx, $rb_src_idx,
            $update_even_points, $idx3_is_even, 
            $depletion_handling, $is_weighting_potential, $only2d)


@code_llvm SolidStateDetectors.innerloop!(line_weights, pssrb, iy, iny, iz, inz, rb_tar_idx, rb_src_idx,
            update_even_points, idx3_is_even, 
            depletion_handling, is_weighting_potential, only2d)

            