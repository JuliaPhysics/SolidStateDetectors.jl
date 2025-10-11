"""
    struct DielectricDistribution{T, N, S, AT} <: AbstractArray{T, N}
        
Dielectric distribution, or distribution of the relative permittivity, ``\\epsilon_r``, 
needed to calculate the [`ElectricPotential`](@ref).
The dielectric distribution is unitless.
        
## Parametric types 
* `T`: Element type of `data`.
* `N`: Dimension of the `grid` and `data` array.  
* `S`: Coordinate system (`Cartesian` or `Cylindrical`).
* `AT`: Axes type.  
        
## Fields
* `data::Array{T, N}`: Array containing the values of the dielectric distribution at the discrete points of the `grid`.
* `grid::Grid{T, N, S, AT}`: [`Grid`](@ref) defining the discrete points at which the electric potential is determined.

!!! note
    The `data` array contains the values of the dielectric distribution at the discrete points 
    between the points defined by the axes ticks of the extended grid of `grid`.\n
    Thus, `size(data) == size(grid) .+ 1` !
"""
struct DielectricDistribution{T, N, S, AT} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S, AT}
end

@inline size(ϵ::DielectricDistribution{T, N, S}) where {T, N, S} = size(ϵ.data)
@inline length(ϵ::DielectricDistribution{T, N, S}) where {T, N, S} = length(ϵ.data)
@inline getindex(ϵ::DielectricDistribution{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(ϵ.data, I...)
@inline getindex(ϵ::DielectricDistribution{T, N, S}, i::Int) where {T, N, S} = getindex(ϵ.data, i)
@inline getindex(ϵ::DielectricDistribution{T, N, S}, s::Symbol) where {T, N, S} = getindex(ϵ.grid, s)



function Base.NamedTuple(ϵ::DielectricDistribution{T}) where {T <: SSDFloat}
    return (
        grid = NamedTuple(ϵ.grid),
        values = ϵ.data,
    )
end
Base.convert(T::Type{NamedTuple}, x::DielectricDistribution) = T(x)
    
function DielectricDistribution(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(ustrip(nt.values[1]))
    S = get_coordinate_system(grid)
    N = get_number_of_dimensions(grid)
    DielectricDistribution{T, N, S, typeof(grid.axes)}( nt.values, grid )
end
Base.convert(T::Type{DielectricDistribution}, x::NamedTuple) = T(x)
