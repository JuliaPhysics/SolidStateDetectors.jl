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

pssrb_car = SolidStateDetectors.PotentialCalculationSetup( 
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


# calculate_weights_for_innerloop! & innerloop! 
depletion_handling = Val{false}()

idx3_is_even = Val(iseven(idx3))

iz = idx3
inz = iz - 1 
            
pwwzr        = pssrb_car.geom_weights[3][1, inz]
pwwzl        = pssrb_car.geom_weights[3][2, inz]
pwΔmpz       = pssrb_car.geom_weights[3][3, inz]
Δz_ext_inv_r = pssrb_car.geom_weights[3][4, inz + 1]
Δz_ext_inv_l = pssrb_car.geom_weights[3][4, inz]

line_weights = Array{T, 2}(undef, size(pssrb_car.potential, 1) - 2, 6)

iy = 2
iny = iy - 1

pwwyr  = pssrb_car.geom_weights[2][1, iny]
pwwyl  = pssrb_car.geom_weights[2][2, iny]
pwΔmpy = pssrb_car.geom_weights[2][3, iny] 
pwΔmpy_pwΔmpz = pwΔmpy * pwΔmpz
Δy_ext_inv_r_pwΔmpz  = pssrb_car.geom_weights[2][4, iny + 1] * pwΔmpz
Δy_ext_inv_l_pwΔmpz  = pssrb_car.geom_weights[2][4, iny]     * pwΔmpz
Δz_ext_inv_r_pwΔmpy = Δz_ext_inv_r * pwΔmpy
Δz_ext_inv_l_pwΔmpy = Δz_ext_inv_l * pwΔmpy

pwwyr_pwwzr = pwwyr * pwwzr
pwwyr_pwwzl = pwwyr * pwwzl
pwwyl_pwwzr = pwwyl * pwwzr
pwwyl_pwwzl = pwwyl * pwwzl


begin
    @info "Grid size: $(size(sim.electric_potential.data))"
    # SolidStateDetectors.outerloop!(pssrb_car, nthreads, update_even_points, depletion_handling, is_weighting_potential, only2d)
    for dp in (false, true)
        depletion_handling = Val{dp}()
        @info "Depletion handling: $dp"
        for nt in (1, 2, 4, 8, 16, 32, 64)
            @info "N Threads: $nt"
            @btime SolidStateDetectors.outerloop!($pssrb_car, $nt, $update_even_points, $depletion_handling, $is_weighting_potential, $only2d)
        end
    end
end


SolidStateDetectors.middleloop!(idx3, rb_tar_idx, rb_src_idx, pssrb_car, 
                        update_even_points, depletion_handling, is_weighting_potential, only2d, idx3iseven)

@code_warntype SolidStateDetectors.middleloop!(idx3, rb_tar_idx, rb_src_idx, pssrb_car, 
                        update_even_points, depletion_handling, is_weighting_potential, only2d, idx3iseven)


@btime SolidStateDetectors.middleloop!($idx3, $rb_tar_idx, $rb_src_idx, $pssrb_car, 
                        $update_even_points, $depletion_handling, $is_weighting_potential, $only2d, $idx3iseven)



SolidStateDetectors.calculate_weights_for_innerloop!(line_weights, pssrb_car, iy, iny, iz, inz,
            update_even_points, idx3_is_even, 
            pwwzr, pwwzl, pwwyr, pwwyl,
            pwwyr_pwwzr, pwwyr_pwwzl, pwwyl_pwwzr, pwwyl_pwwzl,
            pwΔmpy_pwΔmpz,
            Δz_ext_inv_r_pwΔmpy, Δz_ext_inv_l_pwΔmpy,
            Δy_ext_inv_r_pwΔmpz, Δy_ext_inv_l_pwΔmpz) 


@btime SolidStateDetectors.calculate_weights_for_innerloop!($line_weights, $pssrb_car, $iy, $iny, $iz, $inz,
            $update_even_points, $idx3_is_even, 
            $pwwzr, $pwwzl, $pwwyr, $pwwyl,
            $pwwyr_pwwzr, $pwwyr_pwwzl, $pwwyl_pwwzr, $pwwyl_pwwzl,
            $pwΔmpy_pwΔmpz,
            $Δz_ext_inv_r_pwΔmpy, $Δz_ext_inv_l_pwΔmpy,
            $Δy_ext_inv_r_pwΔmpz, $Δy_ext_inv_l_pwΔmpz)

SolidStateDetectors.innerloop!(line_weights, pssrb_car, iy, iny, iz, inz, rb_tar_idx, rb_src_idx,
            update_even_points, idx3_is_even, 
            depletion_handling, is_weighting_potential, only2d)

@code_warntype SolidStateDetectors.innerloop!(line_weights, pssrb_car, iy, iny, iz, inz, rb_tar_idx, rb_src_idx,
            update_even_points, idx3_is_even, 
            depletion_handling, is_weighting_potential, only2d)

@btime SolidStateDetectors.innerloop!($line_weights, $pssrb_car, $iy, $iny, $iz, $inz, $rb_tar_idx, $rb_src_idx,
            $update_even_points, $idx3_is_even, 
            $depletion_handling, $is_weighting_potential, $only2d)


@code_llvm SolidStateDetectors.innerloop!(line_weights, pssrb_car, iy, iny, iz, inz, rb_tar_idx, rb_src_idx,
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

#######################

using KernelAbstractions
using CUDAKernels
using CUDAKernels.CUDA: CuArray

abstract type T32 end
abstract type T64 end

@kernel function test_kernel!(
    a::AbstractArray{T, 2}
) where {T, S}
    linear_idx = @index(Global)
    if linear_idx <= length(a)
        a[linear_idx] = test_dispatch_func(T32)
    end
end

@inline test_dispatch_func(::Type{T32}) = 32f0
@inline test_dispatch_func(::Type{T64}) = 64e0

ca_32 = CuArray(zeros(Float32, 2, 3));
ca_64 = CuArray(zeros(Float64, 2, 3));

k = test_kernel!( CUDAKernels.CUDADevice() )
wait(k(ca_32, ndrange=size(ca_32)))
wait(k(ca_64, ndrange=size(ca_64)))
ca_32[1] == test_dispatch_func(T32)
ca_64[1] == test_dispatch_func(T64)

@kernel function test_kernel!(
    a::AbstractArray{T, 2},
    ::Type{S}
) where {T, S}
    linear_idx = @index(Global)
    if linear_idx <= length(a)
        a[linear_idx] = test_dispatch_func(S)
    end
end

#######################

using Adapt
using KernelAbstractions
using CUDAKernels
using CUDAKernels.CUDA: CuArray
using SolidStateDetectors
using SolidStateDetectors: PotentialCalculationSetup, _guess_optimal_number_of_threads_for_SOR,
    get_device, Cartesian, Cylindrical, PointType, sor_kernel

T = Float32
sim_cpu = Simulation{T}(SSD_examples[:CGD]);
calculate_electric_potential!(sim_cpu, device_array_type = Array, 
    convergence_limit = T(0), refinement_limits = missing,
    max_n_iterations = 2000 )
using Plots; plot(sim_cpu.electric_potential, y = 0)

sim = Simulation{T}(SSD_examples[:CGD]);
calculate_electric_potential!(sim, device_array_type = CuArray, 
    convergence_limit = T(0), refinement_limits = missing,
    max_n_iterations = 2000 )

using Plots; plot(sim.electric_potential, y = 0)

device_array_type = CuArray
CS = SolidStateDetectors.get_coordinate_system(sim)
sor_consts = T(1.2)
not_only_paint_contacts = true
paint_contacts = true
DAT = device_array_type

pssrb = adapt(device_array_type, PotentialCalculationSetup(
    sim.detector, sim.electric_potential.grid, sim.medium, sim.electric_potential.data, sor_consts = T.(sor_consts),
    use_nthreads = _guess_optimal_number_of_threads_for_SOR(size(sim.electric_potential.grid), Base.Threads.nthreads(), CS),    
    not_only_paint_contacts = not_only_paint_contacts, paint_contacts = paint_contacts,
));

update_even_points = true
depletion_handling_enabled = true
_is_weighting_potential = false
only_2d = false

device = get_device(DAT)
N_grid_points = prod(size(pssrb.potential)[1:3] .- 2)

begin
    @kernel function kernelXY001(
        potential::AbstractArray{T, 4},
        point_types::AbstractArray{PointType, 4},
        volume_weights::AbstractArray{T, 4},
        q_eff_imp::AbstractArray{T, 4},
        q_eff_fix::AbstractArray{T, 4},
        ϵ_r::AbstractArray{T, 3},
        geom_weights::NTuple{3, <:AbstractArray{T, 2}},
        sor_const::AbstractArray{T, 1},
        update_even_points::Bool,
        depletion_handling_enabled::Bool,
        is_weighting_potential::Bool,
        only2d::Bool
    ) where {T}
        linear_idx = @index(Global)
        # car_idx = @index(Global, Cartesian)
        # if linear_idx < 20 
        #     @print(linear_idx, " ", car_idx[1], " ", car_idx[2], " ", car_idx[3], "\n")
        # end
        sor_kernel(
            potential,
            point_types,
            volume_weights,
            q_eff_imp,
            q_eff_fix,
            ϵ_r,
            geom_weights,
            sor_const,
            update_even_points,
            depletion_handling_enabled,
            is_weighting_potential,
            only2d, 
            Cartesian,
            linear_idx
        )
    end 
    @info "1"
    kernel = kernelXY001(device)
    @info "2"
    wait(kernel( 
        pssrb.potential, pssrb.point_types, pssrb.volume_weights, pssrb.q_eff_imp, pssrb.q_eff_fix, pssrb.ϵ_r,
        pssrb.geom_weights, pssrb.sor_const, update_even_points, depletion_handling_enabled, _is_weighting_potential, only_2d, 
        # ndrange = size(pssrb.potential)[1:3] .- 2 # Cartesian indexing
        ndrange = N_grid_points   # Linear indexing
    ))
    @info "3"
end