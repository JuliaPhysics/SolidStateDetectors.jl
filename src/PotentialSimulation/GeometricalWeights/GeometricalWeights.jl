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

struct GeometricalAzimutalAxisWeights{T} <: AbstractGeometricalAxisWeights{T}
    weights::Array{T, 2}
end

struct GeometricalRadialAxisWeights{T} <: AbstractGeometricalAxisWeights{T}
    weights::Array{T, 2}
end

function GeometricalCartesianAxisWeights( ax::DiscreteAxis{T, BL, BR} )::GeometricalCartesianAxisWeights{T} where {T, BL, BR}
    axticks::Vector{T} = collect(ax.ticks)
    ax_ext::Vector{T} = get_extended_ticks(ax)
    Δax_ext::Vector{T} = diff(ax_ext)
    Δax_ext_inv::Vector{T} = inv.(Δax_ext)
    mpax::Vector{T} = midpoints(ax_ext)
    Δmpax::Vector{T} = diff(mpax)
    Δmpax_inv::Vector{T} = inv.(Δmpax)
    Δhmpaxr::Vector{T} = mpax[2:end] - axticks
    Δhmpaxl::Vector{T} = axticks - mpax[1:end - 1]
    waxr::Vector{T} = Δmpax_inv .* Δhmpaxr
    waxl::Vector{T} = Δmpax_inv .* Δhmpaxl

    w::Array{T, 2} = zeros(T, 4, length(waxr) + 1)
    w[1, 1:length(waxr)] = waxr
    w[2, 1:length(waxr)] = waxl
    w[3, 1:length(waxr)] = Δmpax
    w[4, :] = Δax_ext_inv

    return GeometricalCartesianAxisWeights{T}( w )
end

