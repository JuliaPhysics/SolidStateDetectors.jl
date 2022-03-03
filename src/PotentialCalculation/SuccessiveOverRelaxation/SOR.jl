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
Eqs. 5.31 to 5.37 in https://mediatum.ub.tum.de/doc/1620954/1620954.pdf.
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
        imp_scale::T,
        neighbor_potentials::NTuple{6,T}, 
        q_eff_imp::T, 
        volume_weight::T,
        sor_const::T
    )::Tuple{T, PointType} where {T}

This function handles the grid points with volumes which are not fully depleted.
The decision depends on the new calculated potential value, `new_potential`, 
for the respective grid point in that iteration and 
the potential values of the neighboring grid points, `neighbor_potentials`.

If `vmin = minimum(neighbor_potentials) < new_potential < vmax = maximum(neighbor_potentials)` => fully depleted

Else => (partially) undepleted => Set `new_potential` to the respective extremum `vext`, either `vmin` or `vmax`.
Depending on how much the potential overshooted `abs(new_potential - vext)`,
the impurity density in the grid point is scaled down via reducing `imp_scale ∈ [0, 1]` of this grid point.

This decision is based on the fact that the potential inside a solid-state 
detector increases monotonically from a `p+`-contact towards an `n+`-contact.
Thus, there cannot be local extrema. 

!!! note
    If a fixed charge impurity is present, e.g. due to a charged passivated surface,
    this handling is probably not valid anymore as the potential between a 
    `p+`-contact towards an `n+`-contact is not required to change monotonically anymore.
    
"""
@inline @fastmath function handle_depletion(
    new_potential::T, 
    imp_scale::T,
    neighbor_potentials::NTuple{6,T}, 
    q_eff_imp::T, 
    volume_weight::T,
    sor_const::T
)::Tuple{T, T} where {T}
    vmin::T = min(neighbor_potentials...)
    vmax::T = max(neighbor_potentials...)
    
    is_ptype = q_eff_imp < 0 
    neighbor_relevant_extremum = is_ptype ? vmin : vmax
    imp_contribution = q_eff_imp * volume_weight * sor_const 
    
    if !iszero(imp_contribution) 
        new_potential_proposal = new_potential - (T(1) - imp_scale) * imp_contribution
        imp_scale -= (new_potential_proposal - neighbor_relevant_extremum) / imp_contribution
        imp_scale = min(max(imp_scale, T(0)), T(1))
        new_potential -= (T(1) - imp_scale) * imp_contribution 
    end

    new_potential, imp_scale
end

include("CPU_outerloop.jl")
include("CPU_middleloop.jl")
include("CPU_innerloop.jl")
include("CPU_update.jl")

include("GPU_kernel_Cylindrical.jl")
include("GPU_kernel_Cartesian3D.jl")
include("GPU_kernel.jl")
include("GPU_update.jl")

include("convergence.jl")


function mark_undep_bits_and_contact_imp_scale!(point_types::Array{PointType, 3}, imp_scale::Array{T, 3}) where {T}
    @inbounds for i in eachindex(imp_scale)
        imp_scale[i] *= !is_fixed_point_type(point_types[i])
        point_types[i] += undepleted_bit * (is_pn_junction_point_type(point_types[i]) && imp_scale[i] < 1)
    end
    nothing
end