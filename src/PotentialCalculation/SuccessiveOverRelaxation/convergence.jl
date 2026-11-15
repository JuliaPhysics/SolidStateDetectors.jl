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
    exit_reason = :converged
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
                exit_reason = :stagnated
                break
            end
        end
        if max_n_iterations != -1 && n_performed_iterations >= max_n_iterations
            exit_reason = :max_iterations
            break
        end
    end
    if exit_reason == :stagnated
        # Two unrelated causes can both land here, so each gets its own _id:
        # maxlog=1 dedupes per _id, not per message content, and without
        # separate ids the first cause encountered in a process would use up
        # the one allowed warning and permanently silence the other.
        if depletion_handling_enabled
            @warn "SOR stopped early: the residual stagnated at c = $c (limit = $c_limit) after $n_performed_iterations iterations without reaching the requested convergence limit. This can happen for undepleted detectors, where grid points near the depletion boundary keep switching classification between iterations." maxlog=1 _id=:sor_stagnated_depletion_handling
        else
            @warn "SOR stopped early: the residual stagnated at c = $c (limit = $c_limit) after $n_performed_iterations iterations without reaching the requested convergence limit. This can happen simply because the requested convergence_limit is close to or below the numerical precision of $T (e.g. eps(Float32) ≈ 1.2e-7)." maxlog=1 _id=:sor_stagnated_precision_limit
        end
    elseif exit_reason == :max_iterations
        @warn "SOR reached max_n_iterations = $max_n_iterations before converging (c = $c, limit = $c_limit). The result may not be fully converged." maxlog=1
    end
    if verbose ProgressMeter.finish!(prog) end
    return c
end
