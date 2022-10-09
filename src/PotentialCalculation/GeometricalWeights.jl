abstract type AbstractGeometricalAxisWeights{T} <: AbstractArray{T, 2} end

function sizeof(gw::AbstractGeometricalAxisWeights{T})::Int where {T}
    return sizeof(gw.weights)
end
@inline function size(gw::AbstractGeometricalAxisWeights{T}) where {T}
    return size(gw.weights)
end

@inline getindex(gw::AbstractGeometricalAxisWeights{T}, I::Vararg{Int, N}) where {T, N} = getindex(gw.weights, I...)

"""
    struct GeometricalCartesianAxisWeights{T} <: AbstractGeometricalAxisWeights{T}

`GeometricalCartesianAxisWeights` stores certain precalculated (in its construction) fixed values
of one of the cartesian (linear) axes, `ax`, of the grid on which the field calculation is performed. 
These precalculated values are further used in the field calculation 
to calculate the six weights for the neighboring grid points in the SOR. 

`GeometricalCartesianAxisWeights` has only one field: `weights::Array{T}, 2`, which is 
of size `(4, length(ax.ticks))`. Each row holding the 4 precalculated values for one axis tick.

Axis ticks, `t = ax.ticks`, and midpoints, `mp = midpoints(get_extended_ticks(ax))`,\\
(the middle points between to axis ticks) inbetween:\\
`... t[i-1] --- mp[i] --- t[i] --- mp[i+1] --- t[i+1] ...`\\
are required for the understanding of the precalculated values.

The columns store the following quantities:
* `weights[1, i]`: `(mp[i+1] - t[i])` / `weights[3, i]`
* `weights[2, i]`: `(t[i] - mp[i])` / `weights[3, i]`
* `weights[3, i]`: `mp[i+1] - mp[i]`
* `weights[4, i]`: `inv(t[i] - t[i-1])`
"""
struct GeometricalCartesianAxisWeights{T} <: AbstractGeometricalAxisWeights{T}
    weights::Array{T, 2}
end


"""
    struct GeometricalAzimutalAxisWeights{T} <: AbstractGeometricalAxisWeights{T}

`GeometricalAzimutalAxisWeights` stores certain precalculated (in its construction) fixed values
of the azimulal (or polar) axis, `ax`, of an cylindrical grid on which the field calculation is performed. 
These precalculated values are further used in the field calculation 
to calculate the six weights for the neighboring grid points in the SOR. 

`GeometricalAzimutalAxisWeights` has only one field: `weights::Array{T}, 2`, which is 
of size `(4, length(ax.ticks))`. Each row holding the 4 precalculated values for one axis tick.

Axis ticks, `t = ax.ticks`, and midpoints, `mp = midpoints(get_extended_ticks(ax))`,\\
(the middle points between to axis ticks) inbetween:\\
`... t[i-1] --- mp[i] --- t[i] --- mp[i+1] --- t[i+1] ...`\\
are required for the understanding of the precalculated values.

The columns store the following quantities:
* `weights[1, i]`: `(mp[i+1] - t[i])` / `weights[3, i]`
* `weights[2, i]`: `(t[i] - mp[i])` / `weights[3, i]`
* `weights[3, i]`: `mp[i+1] - mp[i]`
* `weights[4, i]`: `inv(t[i] - t[i-1])`
* `weights[5, i]`: `normalize(wφ[3, 1:2:end], 1)`
* `weights[6, i]`: `normalize(wφ[3, 2:2:end], 1)`

Developer note: This actually seems to be the same as `GeometricalCartesianAxisWeights`.
Thus, we actually could remove that. Or change this into an alias.
"""
struct GeometricalAzimutalAxisWeights{T} <: AbstractGeometricalAxisWeights{T}
    weights::Array{T, 2}
end

"""
    struct GeometricalRadialAxisWeights{T} <: AbstractGeometricalAxisWeights{T}

`GeometricalRadialAxisWeights` stores certain precalculated (in its construction) fixed values
of the radial axis, `ax`, of an cylindrical grid on which the field calculation is performed. 
These precalculated values are further used in the field calculation 
to calculate the six weights for the neighboring grid points in the SOR. 

`GeometricalRadialAxisWeights` has only one field: `weights::Array{T}, 2`, which is 
of size `(6, length(ax.ticks))`. Each row holding the 6 precalculated values for one axis tick.

Axis ticks, `t = ax.ticks`, and midpoints, `mp = midpoints(get_extended_ticks(ax))`,\\
(the middle points between to axis ticks) inbetween:\\
`... t[i-1] --- mp[i] --- t[i] --- mp[i+1] --- t[i+1] ...`\\
are required for the understanding of the precalculated values.

The columns store the following quantities:
* `weights[1, i]`: `(mp[i+1] - t[i])` / `weights[3, i]`
* `weights[2, i]`: `(t[i] - mp[i])` / `weights[3, i]`
* `weights[3, i]`: `(mp[i+1] - mp[i]) / t[i]`
* `weights[4, i]`: `mp[i+1] / (t[i+1] - t[i])`
* `weights[5, i]`: `mp[i] / (t[i] - t[i-1])`\\
  # Developer note: maybe this one could be removed as it basically stores the same as `weights[4, i-1]`
* `weights[6, i]`: `(mp[i+1]^2 - mp[i]^2)/2`

Developer note: Some of the weights for `i == 1` are set manually as `r = 0` requires some special treatment.
"""
struct GeometricalRadialAxisWeights{T} <: AbstractGeometricalAxisWeights{T}
    weights::Array{T, 2}
end



