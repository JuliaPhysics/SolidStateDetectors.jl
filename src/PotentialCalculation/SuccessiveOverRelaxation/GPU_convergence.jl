function _update_till_convergence!( pssrb::PotentialCalculationSetup{T, S, 3}, 
                                    convergence_limit,
                                    device_array_type::Type{DAT};
                                    n_iterations_between_checks = 500,
                                    depletion_handling::Val{depletion_handling_enabled} = Val{false}(),
                                    only2d::Val{only_2d} = Val{false}(), 
                                    is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
                                    use_nthreads::Int = Base.Threads.nthreads(), 
                                    max_n_iterations::Int = 10_000, # -1
                                    verbose::Bool = true
                                ) where {T, S, depletion_handling_enabled, only_2d, _is_weighting_potential, DAT <: GPUArrays.AbstractGPUArray}
    device = get_device(DAT)
    rb_half_inner_size = size(pssrb.potential)[1:3] .- 2
    ndrange = prod(rb_half_inner_size)
    kernel = get_sor_kernel(S, device)
    c_limit = _is_weighting_potential ? convergence_limit : abs(convergence_limit * pssrb.bias_voltage)
    c = c_limit * 10
    n_performed_iterations = 0
    tmp_potential = similar(pssrb.potential, rb_half_inner_size)
    n_half_grid_points = length(tmp_potential)
    inner_ranges = broadcast(i -> 2:size(tmp_potential, i) + 1, (1, 2, 3))
    is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")
    if verbose prog = ProgressThresh(c; dt = 0.1, desc = "Convergence: ", output = stderr, enabled = !is_logging(stderr)) end
    while c > c_limit
        for k in 1:n_iterations_between_checks-1
            update!(pssrb, kernel, ndrange, depletion_handling_enabled, _is_weighting_potential, only2d)
            n_performed_iterations += 1
        end
        begin
            tmp_potential[:, :, :] .= view(pssrb.potential, inner_ranges..., 1)
            update!(pssrb, kernel, ndrange, depletion_handling_enabled, _is_weighting_potential, only2d)
            tmp_potential[:, :, :] .-= view(pssrb.potential, inner_ranges..., 1)
            c = sum(abs.(tmp_potential)) / n_half_grid_points
            if verbose ProgressMeter.update!(prog, c) end
            n_performed_iterations += 1
        end
        if n_performed_iterations >= max_n_iterations break end
    end
    if verbose ProgressMeter.finish!(prog) end
    return c
end                

@inline function update!(
    pssrb::PotentialCalculationSetup{T, S, 3},
    kernel, ndrange,
    depletion_handling_enabled, 
    _is_weighting_potential, 
    only2d::Val{only_2d}
) where {T, S, only_2d}
    update_even_points = true
    wait(kernel( 
        pssrb.potential, pssrb.point_types, pssrb.volume_weights, pssrb.q_eff_imp, pssrb.q_eff_fix, pssrb.ϵ_r,
        pssrb.geom_weights, pssrb.sor_const, update_even_points, depletion_handling_enabled, _is_weighting_potential, only_2d, 
        ndrange = ndrange
    ))
    apply_boundary_conditions!(pssrb, Val(update_even_points), only2d)
    update_even_points = false
    wait(kernel( 
        pssrb.potential, pssrb.point_types, pssrb.volume_weights, pssrb.q_eff_imp, pssrb.q_eff_fix, pssrb.ϵ_r,
        pssrb.geom_weights, pssrb.sor_const, update_even_points, depletion_handling_enabled, _is_weighting_potential, only_2d,
        ndrange = ndrange
    ))
    apply_boundary_conditions!(pssrb, Val(update_even_points), only2d)
    return nothing
end
