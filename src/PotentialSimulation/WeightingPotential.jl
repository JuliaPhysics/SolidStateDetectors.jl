"""
    struct WeightingPotential{T, N, S, AT} <: AbstractArray{T, N}
        
Weighting potential for a given contact which is a unitless quantity.
        
## Parametric types 
* `T`: Element type of `data`.
* `N`: Dimension of the `grid` and `data` array.  
* `S`: Coordinate system (`Cartesian` or `Cylindrical`).
* `AT`: Axes type.  
        
## Fields
* `data::Array{T, N}`: Array containing the values of the weighting potential at the discrete points of the `grid`.
* `grid::Grid{T, N, S, AT}`: Grid defining the discrete points for which the weighting potential is determined.
"""
struct WeightingPotential{T, N, S <: AbstractCoordinateSystem, AT} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S, AT}
end

@inline size(wp::WeightingPotential{T, N, S}) where {T, N, S} = size(wp.data)
@inline length(wp::WeightingPotential{T, N, S}) where {T, N, S} = length(wp.data)
@inline getindex(wp::WeightingPotential{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(wp.data, I...)
@inline getindex(wp::WeightingPotential{T, N, S}, i::Int) where {T, N, S} = getindex(wp.data, i)
@inline getindex(wp::WeightingPotential{T, N, S}, s::Symbol) where {T, N, S} = getindex(wp.grid, s)


function WeightingPotential(fss::PotentialSimulationSetup{T, 3, Cylindrical}; kwargs...)::WeightingPotential{T, 3, Cylindrical} where {T <: SSDFloat}
    return get_2π_potential(WeightingPotential{T, 3, Cylindrical}(fss.potential, fss.grid); kwargs...)
end
function WeightingPotential(fss::PotentialSimulationSetup{T, 3, Cartesian})::WeightingPotential{T, 3, Cartesian} where {T <: SSDFloat}
    return WeightingPotential{T, 3, Cartesian}(fss.potential, fss.grid)
end

function WeightingPotential(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(ustrip(nt.values[1]))
    S = get_coordinate_system(grid)
    N = get_number_of_dimensions(grid)
    WeightingPotential{T, N, S}( nt.values, grid)
end
Base.convert(T::Type{WeightingPotential}, x::NamedTuple) = T(x)

function NamedTuple(ep::WeightingPotential{T, 3}) where {T}
    return (
        grid = NamedTuple(ep.grid),
        values = ep.data,
    )
end
Base.convert(T::Type{NamedTuple}, x::WeightingPotential) = T(x)
