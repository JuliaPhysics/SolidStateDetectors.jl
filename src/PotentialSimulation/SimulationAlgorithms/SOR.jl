function update!(   pssrb::PotentialSimulationSetupRB{T}; use_nthreads::Int = Base.Threads.nthreads(), 
        depletion_handling::Val{depletion_handling_enabled} = Val{false}(), only2d::Val{only_2d} = Val{false}(),
        is_weighting_potential::Val{_is_weighting_potential} = Val{false}())::Nothing where {T, depletion_handling_enabled, only_2d, _is_weighting_potential}
    outerloop!(pssrb, use_nthreads, Val{true}(), depletion_handling, is_weighting_potential, only2d)
    apply_boundary_conditions!(pssrb, Val{true}(), only2d)
    outerloop!(pssrb, use_nthreads, Val{false}(), depletion_handling, is_weighting_potential, only2d)
    apply_boundary_conditions!(pssrb, Val{false}(), only2d)
    nothing
end

@fastmath function outerloop!( pssrb::PotentialSimulationSetupRB{T}, use_nthreads::Int,
    update_even_points::Val{even_points},
    depletion_handling::Val{depletion_handling_enabled},
    is_weighting_potential::Val{_is_weighting_potential},
    only2d::Val{only_2d})::Nothing where {T, even_points, depletion_handling_enabled, _is_weighting_potential, only_2d}
    @inbounds begin 
        rb_tar_idx::Int, rb_src_idx::Int = even_points ? (rb_even::Int, rb_odd::Int) : (rb_odd::Int, rb_even::Int) 

        @onthreads 1:use_nthreads for idx3 in workpart(2:1:(size(pssrb.potential, 3) - 1), 1:use_nthreads, Base.Threads.threadid())
            middleloop!( idx3, rb_tar_idx, rb_src_idx, pssrb, 
                         update_even_points, depletion_handling, 
                         is_weighting_potential, only2d, Val(iseven(idx3))::Union{Val{true}, Val{false}})
        end 
    end 
    nothing
end


function _guess_optimal_number_of_threads_for_SOR(gs::NTuple{3, Integer}, max_nthreads::Integer, S::Union{Type{Cylindrical}, Type{Cartesian}})::Int
    max_nthreads = min(Base.Threads.nthreads(), max_nthreads)
    n = S == Cylindrical ? gs[2] * gs[3] : gs[1] * gs[2] # Number of grid points to be updated in each iteration of the outer loop
    #=
        Due to threading overhead calculating small grids with too many threads results in a worse performance.
        Thus, we limit them depending on the grid size.
        These ranges may depend on the hardware and are set by educated guesses.
        # nt = if n < 100
        #     4
        # elseif n < 200
        #     8
        # elseif n < 400
        #     16
        # elseif n < 800
        #     32
        # elseif n < 1600
        #     64
        # else
        #     max_nthreads
        # end
    =#
    return min(nextpow(2, max(cld(n+1, 25), 4)), max_nthreads)
end