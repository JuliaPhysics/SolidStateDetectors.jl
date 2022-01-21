"""
    function update!(   
    pssrb::PotentialSimulationSetupRB{T}; 
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
    pssrb::PotentialSimulationSetupRB{T}; 
    use_nthreads::Int = Base.Threads.nthreads(), 
    depletion_handling::Val{depletion_handling_enabled} = Val{false}(), 
    only2d::Val{only_2d} = Val{false}(),
    is_weighting_potential::Val{_is_weighting_potential} = Val{false}()
)::Nothing where {T, depletion_handling_enabled, only_2d, _is_weighting_potential}
    outerloop!(pssrb, use_nthreads, Val{true}(), depletion_handling, is_weighting_potential, only2d)
    apply_boundary_conditions!(pssrb, Val{true}(), only2d)
    outerloop!(pssrb, use_nthreads, Val{false}(), depletion_handling, is_weighting_potential, only2d)
    apply_boundary_conditions!(pssrb, Val{false}(), only2d)
    nothing
end

function update_and_get_max_abs_diff!(  pssrb::PotentialSimulationSetupRB{T, S, DAT, N1, N2},
                                        depletion_handling::Val{depletion_handling_enabled}, 
                                        only2d::Val{only_2d} = Val{false}(), 
                                        is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
                                        use_nthreads::Int = Base.Threads.nthreads()
                                        )::T where {T, S, DAT, N1, N2, depletion_handling_enabled, only_2d, _is_weighting_potential}
    tmp_potential::Array{T, N2} = copy(pssrb.potential)
    if depletion_handling_enabled
        update!(pssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
        slopes::Array{T, N2} = tmp_potential - pssrb.potential
        @inbounds for i in 1:19
            tmp_potential[:] = pssrb.potential[:]
            update!(pssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
            slopes += tmp_potential - pssrb.potential
        end
        @inbounds slopes /= 20
        return maximum(abs.(slopes))
    else
        for i in 1:10
            update!(pssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
        end
        max_diff::T = maximum(abs.(tmp_potential - pssrb.potential))
        return max_diff
    end
end

function _update_till_convergence!( pssrb::PotentialSimulationSetupRB{T, S, 3}, 
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
    cl::T = _is_weighting_potential ? convergence_limit : abs(convergence_limit * pssrb.bias_voltage) # to get relative change in respect to bias voltage
    # To disable automatically ProgressMeters in CI builds:
    is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")
    if verbose prog = ProgressThresh(cl; dt = 0.1, desc = "Convergence: ", output = stderr, enabled = !is_logging(stderr)) end
    while cf > cl
        for i in 1:n_iterations_between_checks
            update!(pssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
        end
        cf = update_and_get_max_abs_diff!(pssrb, depletion_handling, only2d, is_weighting_potential, use_nthreads)
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
        tmp_point_types = pssrb.point_types .& undepleted_bit
        for i in 1:10
            update!(pssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
            @inbounds for i in eachindex(pssrb.point_types)
                if (pssrb.point_types[i] & undepleted_bit == 0) && (tmp_point_types[i] > 0)
                    pssrb.point_types[i] += undepleted_bit
                elseif (pssrb.point_types[i] & undepleted_bit > 0) && (tmp_point_types[i] == 0)
                    tmp_point_types[i] += undepleted_bit
                end
            end
        end
    end

    if verbose ProgressMeter.finish!(prog) end
    return cf
end



# function load_weights_for_innerloop!(
#     line_weights, pssrb::PotentialSimulationSetupRB{T, Nothing, 3, Array{T, 3}},
#         i2, in2, i3, in3, 
#         update_even_points, i23_is_even_t, 
#         pww3r, pww3l, pww2r, pww2l,
#         pww3r_pww2r, pww3l_pww2r, pww3r_pww2l, pww3l_pww2l,
#         pwΔmp2_pwΔmp3,
#         pwΔmp2r, pwΔmp2l,
#         pwΔmp3r, pwΔmp3l, 
# ) where {T}
#     @fastmath @inbounds @simd ivdep for i1 in 2:(size(pssrb.potential, 1) - 1)
#         in1 = nidx(i1, update_even_points, i23_is_even_t)

#         pww1r        = pssrb.geom_weights[3][1, in1]
#         pww1l        = pssrb.geom_weights[3][2, in1]
#         pwΔmp1       = pssrb.geom_weights[3][3, in1]
#         Δ1_ext_inv_l = pssrb.geom_weights[3][4, in1]
#         Δ1_ext_inv_r = pssrb.geom_weights[3][4, in1 + 1]

#         ϵ_lrr = pssrb.ϵ_r[ in1    ,  i2,  i3 ]
#         ϵ_llr = pssrb.ϵ_r[ in1    , in2,  i3 ]
#         ϵ_lrl = pssrb.ϵ_r[ in1    ,  i2, in3 ] 
#         ϵ_lll = pssrb.ϵ_r[ in1    , in2, in3 ] 
#         ϵ_rrr = pssrb.ϵ_r[ in1 + 1,  i2,  i3 ]
#         ϵ_rlr = pssrb.ϵ_r[ in1 + 1, in2,  i3 ]
#         ϵ_rrl = pssrb.ϵ_r[ in1 + 1,  i2, in3 ]
#         ϵ_rll = pssrb.ϵ_r[ in1 + 1, in2, in3 ]

#         pww2r_pww1r = pww2r * pww1r
#         pww2l_pww1r = pww2l * pww1r
#         pww2r_pww1l = pww2r * pww1l
#         pww2l_pww1l = pww2l * pww1l
#         pww3l_pww1r = pww3l * pww1r
#         pww3r_pww1r = pww3r * pww1r
#         pww3l_pww1l = pww3l * pww1l
#         pww3r_pww1l = pww3r * pww1l

#         w1r =        ϵ_rrr * pww3r_pww2r  
#         w1r = muladd(ϵ_rlr,  pww3r_pww2l, w1r)     
#         w1r = muladd(ϵ_rrl,  pww3l_pww2r, w1r)     
#         w1r = muladd(ϵ_rll,  pww3l_pww2l, w1r)
        
#         w1l =        ϵ_lrr * pww3r_pww2r 
#         w1l = muladd(ϵ_llr,  pww3r_pww2l, w1l)    
#         w1l = muladd(ϵ_lrl,  pww3l_pww2r, w1l)    
#         w1l = muladd(ϵ_lll,  pww3l_pww2l, w1l)
        
#         w2r =        ϵ_rrl * pww3l_pww1r 
#         w2r = muladd(ϵ_rrr,  pww3r_pww1r, w2r)  
#         w2r = muladd(ϵ_lrl,  pww3l_pww1l, w2r)    
#         w2r = muladd(ϵ_lrr,  pww3r_pww1l, w2r) 
        
#         w2l =        ϵ_rll * pww3l_pww1r 
#         w2l = muladd(ϵ_rlr,  pww3r_pww1r, w2l)  
#         w2l = muladd(ϵ_lll,  pww3l_pww1l, w2l)    
#         w2l = muladd(ϵ_llr,  pww3r_pww1l, w2l) 
        
#         w3r =        ϵ_rrr * pww2r_pww1r
#         w3r = muladd(ϵ_rlr,  pww2l_pww1r, w3r)   
#         w3r = muladd(ϵ_lrr,  pww2r_pww1l, w3r)    
#         w3r = muladd(ϵ_llr,  pww2l_pww1l, w3r)
        
#         w3l =        ϵ_rrl * pww2r_pww1r
#         w3l = muladd(ϵ_rll,  pww2l_pww1r, w3l)   
#         w3l = muladd(ϵ_lrl,  pww2r_pww1l, w3l)    
#         w3l = muladd(ϵ_lll,  pww2l_pww1l, w3l) 

#         w1l *= Δ1_ext_inv_l * pwΔmp2_pwΔmp3 
#         w1r *= Δ1_ext_inv_r * pwΔmp2_pwΔmp3 
#         w2l *= pwΔmp3l * pwΔmp1
#         w2r *= pwΔmp3r * pwΔmp1
#         w3l *= pwΔmp2l * pwΔmp1 
#         w3r *= pwΔmp2r * pwΔmp1 

#         line_weights[i1-1, 1] = w1l 
#         line_weights[i1-1, 2] = w1r 
#         line_weights[i1-1, 3] = w2l 
#         line_weights[i1-1, 4] = w2r 
#         line_weights[i1-1, 5] = w3l 
#         line_weights[i1-1, 6] = w3r 
#     end
# end