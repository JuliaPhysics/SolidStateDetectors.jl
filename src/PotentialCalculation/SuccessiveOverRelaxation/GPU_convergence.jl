function _update_till_convergence!( pssrb::PotentialCalculationSetup{T, S, 3}, 
                                    convergence_limit::T,
                                    device_array_type::Type{DAT};
                                    n_iterations_between_checks = 500,
                                    depletion_handling::Val{depletion_handling_enabled} = Val{false}(),
                                    only2d::Val{only_2d} = Val{false}(), 
                                    is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
                                    use_nthreads::Int = Base.Threads.nthreads(), 
                                    max_n_iterations::Int = 10_000, # -1
                                    verbose::Bool = true
                                )::T where {T, S, depletion_handling_enabled, only_2d, _is_weighting_potential, DAT <: GPUArrays.AbstractGPUArray}
    device = get_device(DAT)
    N_grid_points = prod(size(pssrb.potential)[1:3] .- 2)
    kernel = get_sor_kernel(S, device)
    @showprogress for i in 1:max_n_iterations
        update_even_points = true
        wait(kernel( 
            pssrb.potential, pssrb.point_types, pssrb.volume_weights, pssrb.q_eff_imp, pssrb.q_eff_fix, pssrb.ϵ_r,
            pssrb.geom_weights, pssrb.sor_const, update_even_points, depletion_handling_enabled, _is_weighting_potential, only_2d, 
            ndrange = N_grid_points
        ))
        apply_boundary_conditions!(pssrb, Val(update_even_points), only2d)
        update_even_points = false
        wait(kernel( 
            pssrb.potential, pssrb.point_types, pssrb.volume_weights, pssrb.q_eff_imp, pssrb.q_eff_fix, pssrb.ϵ_r,
            pssrb.geom_weights, pssrb.sor_const, update_even_points, depletion_handling_enabled, _is_weighting_potential, only_2d,
            ndrange = N_grid_points
        ))
        apply_boundary_conditions!(pssrb, Val(update_even_points), only2d)
    end
    return 0
end                

get_sor_kernel(::Type{Cylindrical}, args...) = sor_cyl_gpu!(args...)
get_sor_kernel(::Type{Cartesian},   args...) = sor_car_gpu!(args...)

function get_device end