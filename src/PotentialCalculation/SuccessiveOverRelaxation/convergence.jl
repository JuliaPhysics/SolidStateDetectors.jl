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
                                ) where {T, S, depletion_handling_enabled, only_2d, _is_weighting_potential, DAT} #  <: GPUArrays.AbstractGPUArray
    device = get_device(DAT)
    rb_half_inner_size = size(pssrb.potential)[1:3] .- 2
    ndrange = prod(rb_half_inner_size)
    kernel = get_sor_kernel(S, device)
    c_limit = _is_weighting_potential ? convergence_limit : abs(convergence_limit * pssrb.bias_voltage)
    c = (one(c_limit) + c_limit) * 10 # Has to be larger than c_limit at the beginning
    n_performed_iterations = 0
    tmp_potential = similar(pssrb.potential, rb_half_inner_size)
    inner_ranges = broadcast(i -> 2:size(tmp_potential, i) + 1, (1, 2, 3))
    is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")
    cs = fill(c, 4) # 4 is chosen by testing
    if verbose prog = ProgressThresh(c_limit; dt = 0.1, desc = "Convergence: ", output = stderr, enabled = !is_logging(stderr)) end
    while c > c_limit
        for _ in 1:n_iterations_between_checks-1
            update!(pssrb, kernel, ndrange; use_nthreads, depletion_handling, is_weighting_potential, only2d)
            n_performed_iterations += 1
        end
        begin
            tmp_potential[:, :, :] .= view(pssrb.potential, inner_ranges..., 1)
            update!(pssrb, kernel, ndrange; use_nthreads, depletion_handling, is_weighting_potential, only2d)
            tmp_potential[:, :, :] .-= view(pssrb.potential, inner_ranges..., 1)
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
        if n_performed_iterations >= max_n_iterations break end
    end
    if verbose ProgressMeter.finish!(prog) end
    return c
end         

# function update_and_get_max_abs_diff!(  pssrb::PotentialCalculationSetup{T, S, DAT, N1, N2},
#                                         depletion_handling::Val{depletion_handling_enabled}, 
#                                         only2d::Val{only_2d} = Val{false}(), 
#                                         is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
#                                         use_nthreads::Int = Base.Threads.nthreads()
#                                         )::T where {T, S, DAT, N1, N2, depletion_handling_enabled, only_2d, _is_weighting_potential}
#     tmp_potential::Array{T, N2} = copy(pssrb.potential)
#     if depletion_handling_enabled
#         update!(pssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
#         slopes::Array{T, N2} = tmp_potential - pssrb.potential
#         @inbounds for i in 1:19
#             tmp_potential[:] = pssrb.potential[:]
#             update!(pssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
#             slopes += tmp_potential - pssrb.potential
#         end
#         @inbounds slopes /= 20
#         return maximum(abs.(slopes))
#     else
#         for i in 1:10
#             update!(pssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
#         end
#         max_diff::T = maximum(abs.(tmp_potential - pssrb.potential))
#         return max_diff
#     end
# end

# function _update_till_convergence!( pssrb::PotentialCalculationSetup{T, S, 3}, 
#                                     convergence_limit::T, 
#                                     ::Type{Array};
#                                     n_iterations_between_checks = 500,
#                                     depletion_handling::Val{depletion_handling_enabled} = Val{false}(),
#                                     only2d::Val{only_2d} = Val{false}(), 
#                                     is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
#                                     use_nthreads::Int = Base.Threads.nthreads(), 
#                                     max_n_iterations::Int = -1,
#                                     verbose::Bool = true
#                                     )::T where {T, S, depletion_handling_enabled, only_2d, _is_weighting_potential}
#     n_iterations::Int = 0
#     cf::T = Inf
#     cfs::Vector{T} = fill(cf, 10)
#     cl::T = _is_weighting_potential ? convergence_limit : abs(convergence_limit * pssrb.bias_voltage) # to get relative change in respect to bias voltage
#     # To disable automatically ProgressMeters in CI builds:
#     is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")
#     if verbose prog = ProgressThresh(cl; dt = 0.1, desc = "Convergence: ", output = stderr, enabled = !is_logging(stderr)) end
#     while cf > cl
#         for i in 1:n_iterations_between_checks
#             update!(pssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
#         end
#         cf = update_and_get_max_abs_diff!(pssrb, depletion_handling, only2d, is_weighting_potential, use_nthreads)
#         @inbounds cfs[1:end-1] = cfs[2:end]
#         @inbounds cfs[end] = cf
#         slope::T = abs(mean(diff(cfs)))
#         if verbose ProgressMeter.update!(prog, cf) end
#         n_iterations += n_iterations_between_checks
#         if slope < cl
#             # @info "Slope is basically 0 -> Converged: $slope"
#             cf = slope
#         end
#         if max_n_iterations > -1 && n_iterations > max_n_iterations 
#             # @show n_iterations_between_checks
#             verbose && @info "Maximum number of iterations reached. (`n_iterations = $(n_iterations)`)"
#             break
#         end
#     end
#     if depletion_handling_enabled
#         tmp_point_types = pssrb.point_types .& undepleted_bit
#         for i in 1:10
#             update!(pssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
#             @inbounds for i in eachindex(pssrb.point_types)
#                 if (pssrb.point_types[i] & undepleted_bit == 0) && (tmp_point_types[i] > 0)
#                     pssrb.point_types[i] += undepleted_bit
#                 elseif (pssrb.point_types[i] & undepleted_bit > 0) && (tmp_point_types[i] == 0)
#                     tmp_point_types[i] += undepleted_bit
#                 end
#             end
#         end
#     end

#     if verbose ProgressMeter.finish!(prog) end
#     return cf
# end