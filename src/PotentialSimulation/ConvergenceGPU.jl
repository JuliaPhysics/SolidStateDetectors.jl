using .CUDA
using GPUArrays
using CUDAKernels

function _update_till_convergence!( pssrb::PotentialSimulationSetupRB{T, S, 3, DAT}, 
                                    convergence_limit::T;
                                    n_iterations_between_checks = 500,
                                    depletion_handling::Val{depletion_handling_enabled} = Val{false}(),
                                    only2d::Val{only_2d} = Val{false}(), 
                                    is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
                                    use_nthreads::Int = Base.Threads.nthreads(), 
                                    max_n_iterations::Int = 20_000, # -1
                                    verbose::Bool = true
                                    )::T where {T, S, DAT<:GPUArrays.AnyGPUArray, depletion_handling_enabled, only_2d, _is_weighting_potential}
    device = CUDAKernels.CUDADevice() # device = KernelAbstractions.get_device(a)
    N_grid_points = prod(size(pssrb.potential)[1:3] .- 2)
    kernel = sor_gpu!(device, 512, N_grid_points)

    for i in 1:max_n_iterations
        update_even_points = true
        event = kernel( pssrb.potential, pssrb.point_types, pssrb.volume_weights, pssrb.q_eff_imp, pssrb.q_eff_fix, pssrb.ϵ_r,
                        pssrb.geom_weights, pssrb.sor_const, update_even_points, depletion_handling, is_weighting_potential, only_2d)
        wait(event)
        apply_boundary_conditions!(pssrb, Val(update_even_points), only2d)
        update_even_points = false
        event = kernel( pssrb.potential, pssrb.point_types, pssrb.volume_weights, pssrb.q_eff_imp, pssrb.q_eff_fix, pssrb.ϵ_r,
                        pssrb.geom_weights, pssrb.sor_const, update_even_points, depletion_handling, is_weighting_potential, only_2d)
        wait(event)
        apply_boundary_conditions!(pssrb, Val(update_even_points), only2d)
    end

    return 0
end                             
