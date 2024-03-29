using BenchmarkTools
using SolidStateDetectors
using Plots
using Unitful

using SolidStateDetectors: rb_even, rb_odd

### Inner loop of SOR 
T = Float32
sim_car = Simulation{T}(SSD_examples[:CGD]);

calculate_electric_potential!( 
    sim_car, depletion_handling = false,
    convergence_limit = 1e-4,
    max_tick_distance = 1u"mm",
    refinement_limits = [0.2, 0.1],
    # refinement_limits = missing, 
); plot(sim_car.electric_potential, y = 0)

# plot(
#     plot(sim_car.point_types, y = 0),
#     plot(sim_car.electric_potential, y = 0),
# )

pcs_car = SolidStateDetectors.PotentialCalculationSetup( 
    sim_car.detector, 
    sim_car.electric_potential.grid, 
    # Grid(sim_car, max_tick_distance = 1u"mm"), 
    SolidStateDetectors.material_properties[SolidStateDetectors.materials["vacuum"]],
    sim_car.electric_potential.data
);


nthreads = 1 #Base.Threads.nthreads()
is_weighting_potential = Val{false}()
depletion_handling = Val{false}()
only2d = Val{false}()
update_even_points = Val{true}()
update_uneven_points = Val{false}()

idx3 = 2
even_points = true
update_even_points = Val(even_points)
is_weighting_potential = Val(false)
rb_tar_idx, rb_src_idx = even_points ? (rb_even, rb_odd) : (rb_odd, rb_even)
depletion_handling = Val{false}()
idx3iseven = Val(iseven(idx3))
idx3istwo = Val(idx3 == 2)


# calculate_weights_for_innerloop! & innerloop! 
depletion_handling = Val{false}()

idx3_is_even = Val(iseven(idx3))

iz = idx3
inz = iz - 1 
            
pwwzr        = pcs_car.geom_weights[3][1, inz]
pwwzl        = pcs_car.geom_weights[3][2, inz]
pwΔmpz       = pcs_car.geom_weights[3][3, inz]
Δz_ext_inv_r = pcs_car.geom_weights[3][4, inz + 1]
Δz_ext_inv_l = pcs_car.geom_weights[3][4, inz]

line_weights = Array{T, 2}(undef, size(pcs_car.potential, 1) - 2, 6)

iy = 2
iny = iy - 1

pwwyr  = pcs_car.geom_weights[2][1, iny]
pwwyl  = pcs_car.geom_weights[2][2, iny]
pwΔmpy = pcs_car.geom_weights[2][3, iny] 
pwΔmpy_pwΔmpz = pwΔmpy * pwΔmpz
Δy_ext_inv_r_pwΔmpz  = pcs_car.geom_weights[2][4, iny + 1] * pwΔmpz
Δy_ext_inv_l_pwΔmpz  = pcs_car.geom_weights[2][4, iny]     * pwΔmpz
Δz_ext_inv_r_pwΔmpy = Δz_ext_inv_r * pwΔmpy
Δz_ext_inv_l_pwΔmpy = Δz_ext_inv_l * pwΔmpy

pwwyr_pwwzr = pwwyr * pwwzr
pwwyr_pwwzl = pwwyr * pwwzl
pwwyl_pwwzr = pwwyl * pwwzr
pwwyl_pwwzl = pwwyl * pwwzl



SolidStateDetectors.middleloop!(idx3, rb_tar_idx, rb_src_idx, pcs_car, 
                        update_even_points, depletion_handling, is_weighting_potential, only2d, idx3iseven, idx3istwo)

@code_warntype SolidStateDetectors.middleloop!(idx3, rb_tar_idx, rb_src_idx, pcs_car, 
                        update_even_points, depletion_handling, is_weighting_potential, only2d, idx3iseven, idx3istwo)


@btime SolidStateDetectors.middleloop!($idx3, $rb_tar_idx, $rb_src_idx, $pcs_car, 
                        $update_even_points, $depletion_handling, $is_weighting_potential, $only2d, $idx3iseven, $idx3istwo)



SolidStateDetectors.calculate_weights_for_innerloop!(line_weights, pcs_car, iy, iny, iz, inz,
            update_even_points, idx3_is_even, 
            pwwzr, pwwzl, pwwyr, pwwyl,
            pwwyr_pwwzr, pwwyr_pwwzl, pwwyl_pwwzr, pwwyl_pwwzl,
            pwΔmpy_pwΔmpz,
            Δz_ext_inv_r_pwΔmpy, Δz_ext_inv_l_pwΔmpy,
            Δy_ext_inv_r_pwΔmpz, Δy_ext_inv_l_pwΔmpz) 


@btime SolidStateDetectors.calculate_weights_for_innerloop!($line_weights, $pcs_car, $iy, $iny, $iz, $inz,
            $update_even_points, $idx3_is_even, 
            $pwwzr, $pwwzl, $pwwyr, $pwwyl,
            $pwwyr_pwwzr, $pwwyr_pwwzl, $pwwyl_pwwzr, $pwwyl_pwwzl,
            $pwΔmpy_pwΔmpz,
            $Δz_ext_inv_r_pwΔmpy, $Δz_ext_inv_l_pwΔmpy,
            $Δy_ext_inv_r_pwΔmpz, $Δy_ext_inv_l_pwΔmpz)

SolidStateDetectors.innerloop!(line_weights, pcs_car, iy, iny, iz, inz, rb_tar_idx, rb_src_idx,
            update_even_points, idx3_is_even, 
            depletion_handling, is_weighting_potential, only2d)

@code_warntype SolidStateDetectors.innerloop!(line_weights, pcs_car, iy, iny, iz, inz, rb_tar_idx, rb_src_idx,
            update_even_points, idx3_is_even, 
            depletion_handling, is_weighting_potential, only2d)

@btime SolidStateDetectors.innerloop!($line_weights, $pcs_car, $iy, $iny, $iz, $inz, $rb_tar_idx, $rb_src_idx,
            $update_even_points, $idx3_is_even, 
            $depletion_handling, $is_weighting_potential, $only2d)


@code_llvm SolidStateDetectors.innerloop!(line_weights, pcs_car, iy, iny, iz, inz, rb_tar_idx, rb_src_idx,
            update_even_points, idx3_is_even, 
            depletion_handling, is_weighting_potential, only2d)

            
using KernelAbstractions
using CUDAKernels, SolidStateDetectors, Unitful
using CUDAKernels.CUDA: CuArray

T = Float32
sim_car = Simulation{T}(SSD_examples[:CGD]);

calculate_electric_potential!( 
    sim_car, depletion_handling = false,
    device_array_type = CuArray,
    convergence_limit = 1e-4,
    max_tick_distance = 1u"mm",
    refinement_limits = [0.2, 0.1],
    max_n_iterations = 2000,
    # refinement_limits = missing, 
)

begin
    @info "Grid size: $(size(sim.electric_potential.data))"
    # SolidStateDetectors.outerloop!(pcs_car, nthreads, update_even_points, depletion_handling, is_weighting_potential, only2d)
    for dp in (false, true)
        depletion_handling = Val{dp}()
        @info "Depletion handling: $dp"
        for nt in (1, 2, 4, 8, 16, 32, 64)
            @info "N Threads: $nt"
            @btime SolidStateDetectors.outerloop!($pcs_car, $nt, $update_even_points, $depletion_handling, $is_weighting_potential, $only2d)
        end
    end
end
