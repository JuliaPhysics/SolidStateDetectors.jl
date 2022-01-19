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