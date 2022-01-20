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

"""
    function handle_depletion(
        new_potential::T, 
        point_type::PointType, 
        neighbor_potentials::NTuple{6,T}, 
        q_eff_imp::T, 
        volume_weight::T,
        sor_const::T
    )::Tuple{T, PointType} where {T}

This function checks whether the current grid point is depleted or undepleted.
The decision depends on the new calculated potential value, `new_potential`, 
for the respective grid point in that iteration and 
the potential values of the neighboring grid points, `neighbor_potentials`.

If `minimum(neighbor_potentials) < new_potential < maximum(neighbor_potentials)` => depleted

Else => undepleted

This decision is based on the fact that the potential inside a solid-state 
detector increases monotonically from a `p+`-contact towards an `n+`-contact.

!!! note
    If a fixed charge impurity is present, e.g. due to a charged passivated surface,
    this handling is probably not valid anymore as the potential between a 
    `p+`-contact towards an `n+`-contact is not required to change monotonically anymore.
    
"""
@inline @fastmath function handle_depletion(
    new_potential::T, 
    point_type::PointType, 
    neighbor_potentials::NTuple{6,T}, 
    q_eff_imp::T, 
    volume_weight::T,
    sor_const::T
)::Tuple{T, PointType} where {T}
    vmin::T = min(neighbor_potentials...)
    vmax::T = max(neighbor_potentials...)

    if new_potential < vmin || new_potential > vmax
        new_potential -= q_eff_imp * volume_weight * sor_const
        if (point_type & undepleted_bit == 0 && point_type & pn_junction_bit > 0) 
            point_type += undepleted_bit
        end
    elseif point_type & undepleted_bit > 0
        point_type -= undepleted_bit
    end
    new_potential, point_type
end



function _guess_optimal_number_of_threads_for_SOR(gs::NTuple{3, Integer}, max_nthreads::Integer, S::Union{Type{Cylindrical}, Type{Cartesian}})::Int
    max_nthreads = min(Base.Threads.nthreads(), max_nthreads)
    n = S == Cylindrical ? gs[2] * gs[3] : gs[1] * gs[2] # Number of grid points to be updated in each iteration of the outer loop
    return min(nextpow(2, max(cld(n+1, 25), 4)), max_nthreads)
end

include("CPU_outerloop.jl")
include("CPU_middleloop_Cylindrical.jl")
include("CPU_middleloop_Cartesian3D.jl")
include("CPU_innerloop.jl")

include("GPU_Cylindrical.jl")
include("GPU_Cartesian3D.jl")

function update!(   pssrb::PotentialSimulationSetupRB{T}; use_nthreads::Int = Base.Threads.nthreads(), 
    depletion_handling::Val{depletion_handling_enabled} = Val{false}(), only2d::Val{only_2d} = Val{false}(),
    is_weighting_potential::Val{_is_weighting_potential} = Val{false}())::Nothing where {T, depletion_handling_enabled, only_2d, _is_weighting_potential}
outerloop!(pssrb, use_nthreads, Val{true}(), depletion_handling, is_weighting_potential, only2d)
apply_boundary_conditions!(pssrb, Val{true}(), only2d)
outerloop!(pssrb, use_nthreads, Val{false}(), depletion_handling, is_weighting_potential, only2d)
apply_boundary_conditions!(pssrb, Val{false}(), only2d)
nothing
end