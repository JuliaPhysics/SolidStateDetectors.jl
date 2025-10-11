function _update_till_convergence!( pcs::PotentialCalculationSetup{T, S, 3}, 
                                    convergence_limit,
                                    via_KernelAbstractions::Bool;
                                    n_iterations_between_checks = 500,
                                    depletion_handling::Val{depletion_handling_enabled} = Val{false}(),
                                    only2d::Val{only_2d} = Val{false}(), 
                                    is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
                                    use_nthreads::Int = Base.Threads.nthreads(), 
                                    max_n_iterations::Int = 10_000, # -1
                                    verbose::Bool = true
                                ) where {T, S, depletion_handling_enabled, only_2d, _is_weighting_potential} #  <: GPUArrays.AbstractGPUArray
    backend = _ka_get_backend(pcs.potential)
    ndrange = size(pcs.potential)[1:3] .- 2
    kernel = get_sor_kernel(S, backend, Val(via_KernelAbstractions))
    c_limit = _is_weighting_potential ? convergence_limit : abs(convergence_limit * (iszero(pcs.bias_voltage) ? maximum(abs.(pcs.potential)) : pcs.bias_voltage))
    c = (one(c_limit) + c_limit) * 10 # Has to be larger than c_limit at the beginning
    n_performed_iterations = 0
    tmp_potential = similar(pcs.potential, ndrange)
    inner_ranges = broadcast(i -> 2:size(tmp_potential, i) + 1, (1, 2, 3))
    is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")
    cs = fill(c, 4) # 4 is chosen by testing
    if verbose prog = ProgressThresh(c_limit; dt = 0.1, desc = "Convergence: ", output = stderr, enabled = !is_logging(stderr)) end
    while c > c_limit
        for _ in 1:n_iterations_between_checks-1
            update!(pcs, kernel, ndrange; use_nthreads, depletion_handling, is_weighting_potential, only2d)
            n_performed_iterations += 1
        end
        begin
            tmp_potential[:, :, :] .= view(pcs.potential, inner_ranges..., 1)
            update!(pcs, kernel, ndrange; use_nthreads, depletion_handling, is_weighting_potential, only2d)
            tmp_potential[:, :, :] .-= view(pcs.potential, inner_ranges..., 1)
            n_performed_iterations += 1
            c = maximum(abs.(tmp_potential)) 
            if verbose ProgressMeter.update!(prog, c) end
            cs = circshift(cs, -1)
            cs[end] = c
            cs_μ = mean(cs)
            cs_σ = std(cs, mean = cs_μ)
            if cs_σ < c_limit
                # Convergence limit not reached but the value of c does not change anymore
                # Especially needed in case of undepleted detectors (grid points switching between depleted and undepleted)
                break
            end
        end
        if max_n_iterations != -1 && n_performed_iterations >= max_n_iterations break end
    end
    if verbose ProgressMeter.finish!(prog) end
    return c
end         
