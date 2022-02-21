"""
    function update!(   
    pcs::PotentialCalculationSetup{T}; 
    use_nthreads::Int = Base.Threads.nthreads(), 
    depletion_handling::Val{depletion_handling_enabled} = Val{false}(), 
    only2d::Val{only_2d} = Val{false}(),
    is_weighting_potential::Val{_is_weighting_potential} = Val{false}()
)::Nothing where {T, depletion_handling_enabled, only_2d, _is_weighting_potential}

This function performs one iteration of the SOR. One iteration consists out of 4 steps:

    1) Iterate in parallel over all even points and update their potential. 
    2) Apply the boundary conditions at the ends of the grid for all even points. 
    3) Iterate in parallel over all odd points and update their potential. 
    2) Apply the boundary conditions at the ends of the grid for all odd points. 
"""
function update!(   
    pcs::PotentialCalculationSetup{T}; 
    use_nthreads::Int = Base.Threads.nthreads(), 
    depletion_handling::Val{depletion_handling_enabled} = Val{false}(), 
    only2d::Val{only_2d} = Val{false}(),
    is_weighting_potential::Val{_is_weighting_potential} = Val{false}()
)::Nothing where {T, depletion_handling_enabled, only_2d, _is_weighting_potential}
    outerloop!(pcs, use_nthreads, Val{true}(), depletion_handling, is_weighting_potential, only2d)
    apply_boundary_conditions!(pcs, Val{true}(), only2d)
    outerloop!(pcs, use_nthreads, Val{false}(), depletion_handling, is_weighting_potential, only2d)
    apply_boundary_conditions!(pcs, Val{false}(), only2d)
    nothing
end

function update_and_get_max_abs_diff!(  pcs::PotentialCalculationSetup{T, S, DAT, N1, N2},
                                        depletion_handling::Val{depletion_handling_enabled}, 
                                        only2d::Val{only_2d} = Val{false}(), 
                                        is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
                                        use_nthreads::Int = Base.Threads.nthreads()
                                        )::T where {T, S, DAT, N1, N2, depletion_handling_enabled, only_2d, _is_weighting_potential}
    tmp_potential::Array{T, N2} = copy(pcs.potential)
    if depletion_handling_enabled
        update!(pcs, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
        slopes::Array{T, N2} = tmp_potential - pcs.potential
        @inbounds for i in 1:19
            tmp_potential[:] = pcs.potential[:]
            update!(pcs, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
            slopes += tmp_potential - pcs.potential
        end
        @inbounds slopes /= 20
        return maximum(abs.(slopes))
    else
        for i in 1:10
            update!(pcs, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
        end
        max_diff::T = maximum(abs.(tmp_potential - pcs.potential))
        return max_diff
    end
end

function _update_till_convergence!( pcs::PotentialCalculationSetup{T, S, 3}, 
                                    convergence_limit::T, 
                                    ::Type{Array};
                                    n_iterations_between_checks = 500,
                                    depletion_handling::Val{depletion_handling_enabled} = Val{false}(),
                                    only2d::Val{only_2d} = Val{false}(), 
                                    is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
                                    use_nthreads::Int = Base.Threads.nthreads(), 
                                    max_n_iterations::Int = -1,
                                    verbose::Bool = true
                                    )::T where {T, S, depletion_handling_enabled, only_2d, _is_weighting_potential}
    n_iterations::Int = 0
    cf::T = Inf
    cfs::Vector{T} = fill(cf, 10)
    cl::T = _is_weighting_potential ? convergence_limit : abs(convergence_limit * pcs.bias_voltage) # to get relative change in respect to bias voltage
    # To disable automatically ProgressMeters in CI builds:
    is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")
    if verbose prog = ProgressThresh(cl; dt = 0.1, desc = "Convergence: ", output = stderr, enabled = !is_logging(stderr)) end
    while cf > cl
        for i in 1:n_iterations_between_checks
            update!(pcs, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
        end
        cf = update_and_get_max_abs_diff!(pcs, depletion_handling, only2d, is_weighting_potential, use_nthreads)
        @inbounds cfs[1:end-1] = cfs[2:end]
        @inbounds cfs[end] = cf
        slope::T = abs(mean(diff(cfs)))
        if verbose ProgressMeter.update!(prog, cf) end
        n_iterations += n_iterations_between_checks
        if slope < cl
            # @info "Slope is basically 0 -> Converged: $slope"
            cf = slope
        end
        if max_n_iterations > -1 && n_iterations > max_n_iterations 
            # @show n_iterations_between_checks
            verbose && @info "Maximum number of iterations reached. (`n_iterations = $(n_iterations)`)"
            break
        end
    end
    if depletion_handling_enabled
        tmp_point_types = pcs.point_types .& undepleted_bit
        for i in 1:10
            update!(pcs, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
            @inbounds for i in eachindex(pcs.point_types)
                if (pcs.point_types[i] & undepleted_bit == 0) && (tmp_point_types[i] > 0)
                    pcs.point_types[i] += undepleted_bit
                elseif (pcs.point_types[i] & undepleted_bit > 0) && (tmp_point_types[i] == 0)
                    tmp_point_types[i] += undepleted_bit
                end
            end
        end
    end

    if verbose ProgressMeter.finish!(prog) end
    return cf
end