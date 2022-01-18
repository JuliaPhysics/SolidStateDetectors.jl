"""
    function calc_new_potential_SOR_3D(
        q_eff::T,
        volume_weight::T,
        weights::NTuple{6,T},
        neighbor_potentials::NTuple{6,T},
        old_potential::T,
        sor_const::T
    ) where {T}

Calculates and returns the new potential value of a grid point given
its current (old) potential `old_potential`,
the potential values of the six neighbouring grid points (2 in each dimension) `neighbor_potentials`,
the weights corresponding to those six potentials of the neighboring grid points,
the weight due to the volume of the grid point itself (the voxel around it)
and the set SOR constant `sor_const`.

For more detailed information see Chapter 5.2.2 "Calculation of the Electric Potential"
Eqs. 5.36 & 5.37 in https://mediatum.ub.tum.de/doc/1620954/1620954.pdf.
"""
@inline function calc_new_potential_SOR_3D(
    q_eff::T,
    volume_weight::T,
    weights::NTuple{6,T},
    neighbor_potentials::NTuple{6,T},
    old_potential::T,
    sor_const::T
) where {T}
    new_potential = q_eff
    new_potential = muladd(weights[1], neighbor_potentials[1], new_potential)
    new_potential = muladd(weights[2], neighbor_potentials[2], new_potential)
    new_potential = muladd(weights[3], neighbor_potentials[3], new_potential)
    new_potential = muladd(weights[4], neighbor_potentials[4], new_potential)
    new_potential = muladd(weights[5], neighbor_potentials[5], new_potential)
    new_potential = muladd(weights[6], neighbor_potentials[6], new_potential)
    new_potential *= volume_weight
    new_potential -= old_potential
    new_potential = muladd(new_potential, sor_const, old_potential)
    return new_potential
end

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