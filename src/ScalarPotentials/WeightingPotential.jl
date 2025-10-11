"""
    struct WeightingPotential{T, N, S, AT} <: AbstractArray{T, N}
        
Weighting potential for a given [`Contact`](@ref) which is a unitless potential.
        
## Parametric types 
* `T`: Element type of `data`.
* `N`: Dimension of the `grid` and `data` array.  
* `S`: Coordinate system (`Cartesian` or `Cylindrical`).
* `AT`: Axes type.  
        
## Fields
* `data::Array{T, N}`: Array containing the values of the weighting potential at the discrete points of the `grid`.
* `grid::Grid{T, N, S, AT}`: [`Grid`](@ref) defining the discrete points for which the weighting potential is determined.
"""
struct WeightingPotential{T, N, S <: AbstractCoordinateSystem, AT} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S, AT}
end

@inline size(wpot::WeightingPotential{T, N, S}) where {T, N, S} = size(wpot.data)
@inline length(wpot::WeightingPotential{T, N, S}) where {T, N, S} = length(wpot.data)
@inline getindex(wpot::WeightingPotential{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(wpot.data, I...)
@inline getindex(wpot::WeightingPotential{T, N, S}, i::Int) where {T, N, S} = getindex(wpot.data, i)
@inline getindex(wpot::WeightingPotential{T, N, S}, s::Symbol) where {T, N, S} = getindex(wpot.grid, s)


function WeightingPotential(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(ustrip(nt.values[1]))
    S = get_coordinate_system(grid)
    N = get_number_of_dimensions(grid)
    WeightingPotential{T, N, S, typeof(grid.axes)}( nt.values, grid)
end
Base.convert(T::Type{WeightingPotential}, x::NamedTuple) = T(x)

function Base.NamedTuple(wpot::WeightingPotential{T, 3}) where {T}
    return (
        grid = NamedTuple(wpot.grid),
        values = wpot.data,
    )
end
Base.convert(T::Type{NamedTuple}, x::WeightingPotential) = T(x)
