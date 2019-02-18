abstract type AbstractGeometricalAxisWeights{T} <: AbstractArray{T, 2} end

function sizeof(gw::AbstractGeometricalAxisWeights{T})::Int where {T}
    return sizeof(gw.weights)
end
@inline function size(gw::AbstractGeometricalAxisWeights{T}) where {T}
    return size(gw.weights)
end

@inline getindex(gw::AbstractGeometricalAxisWeights{T}, I::Vararg{Int, N}) where {T, N} = getindex(gw.weights, I...)


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

