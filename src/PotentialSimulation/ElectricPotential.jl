"""
    struct ElectricPotential{T, N, S, AT} <: AbstractArray{T, N}
        
Electric potential of the simulation in units of volt (V).
        
## Parametric types 
* `T`: Element type of `data`.
* `N`: Dimension of the `grid` and `data` array.  
* `S`: Coordinate system (`Cartesian` or `Cylindrical`).
* `AT`: Axes type.  
        
## Fields
* `data::Array{T, N}`: Array containing the values of the electric potential at the discrete points of the `grid`.
* `grid::Grid{T, N, S, AT}`: Grid defining the discrete points for which the electric potential is determined.
"""
struct ElectricPotential{T, N, S, AT} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S, AT}
end

@inline size(ep::ElectricPotential{T, N, S}) where {T, N, S} = size(ep.data)
@inline length(ep::ElectricPotential{T, N, S}) where {T, N, S} = length(ep.data)
@inline getindex(ep::ElectricPotential{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(ep.data, I...)
@inline getindex(ep::ElectricPotential{T, N, S}, i::Int) where {T, N, S} = getindex(ep.data, i)
@inline getindex(ep::ElectricPotential{T, N, S}, s::Symbol) where {T, N, S} = getindex(ep.grid, s)


# """
#     ElectricPotential(setup::PotentialSimulationSetup{T, 3, Cylindrical} ; kwargs...)::ElectricPotential{T, 3, Cylindrical}
# 
# Extracts the electric potential from `setup` and extrapolate it to an 2π grid.
# 
# For 2D grids (r and z) the user has to set the keyword `n_points_in_φ::Int`, e.g.: `n_points_in_φ = 36`.
# """
function ElectricPotential(setup::PotentialSimulationSetup{T, 3, Cylindrical} ; kwargs...)::ElectricPotential{T, 3, Cylindrical} where {T}
    return get_2π_potential(ElectricPotential{T, 3, Cylindrical}(setup.potential, setup.grid); kwargs...)
end

# """
#     ElectricPotential(setup::PotentialSimulationSetup{T, 3, Cartesian} ; kwargs...)::ElectricPotential{T, 3, Cartesian}
# 
# Extracts the electric potential from `setup`.
# """
function ElectricPotential(setup::PotentialSimulationSetup{T, 3, Cartesian} )::ElectricPotential{T, 3, Cartesian} where {T}
    return ElectricPotential{T, 3, Cartesian}(setup.potential, setup.grid)
end



function NamedTuple(ep::ElectricPotential{T, 3}) where {T}
    return (
        grid = NamedTuple(ep.grid),
        values = ep.data * u"V",
    )
end
Base.convert(T::Type{NamedTuple}, x::ElectricPotential) = T(x)
    
function ElectricPotential(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(ustrip(nt.values[1]))
    S = get_coordinate_system(grid)
    N = get_number_of_dimensions(grid)
    ElectricPotential{T, N, S}( ustrip.(uconvert.(u"V", nt.values)), grid)
end
Base.convert(T::Type{ElectricPotential}, x::NamedTuple) = T(x)