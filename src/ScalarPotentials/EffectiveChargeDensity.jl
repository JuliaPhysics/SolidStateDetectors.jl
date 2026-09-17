"""
    struct EffectiveChargeDensity{T, N, S, AT} <: AbstractArray{T, N}
        
Effective charge density needed to calculate the [`ElectricPotential`](@ref).
The effective charge density is the charge density (in C/m³) multiplied with the vacuum permittivity 
ϵ0 = 8.854 × 10⁻¹² F/m. Thus, the unit of the effective charge density is Coulomb (V/m²).
        
## Parametric types 
* `T`: Element type of `data`.
* `N`: Dimension of the `grid` and `data` array.  
* `S`: Coordinate system (`Cartesian` or `Cylindrical`).
* `AT`: Axes type.  
        
## Fields
* `data::Array{T, N}`: Array containing the values of the effective charge density at the discrete points of the `grid`.
* `grid::Grid{T, N, S, AT}`: [`Grid`](@ref) defining the discrete points at which the electric potential is determined.
"""
struct EffectiveChargeDensity{T, N, S, AT} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S, AT}
end

@inline size(ρ::EffectiveChargeDensity{T, N, S}) where {T, N, S} = size(ρ.data)
@inline length(ρ::EffectiveChargeDensity{T, N, S}) where {T, N, S} = length(ρ.data)
@inline getindex(ρ::EffectiveChargeDensity{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(ρ.data, I...)
@inline getindex(ρ::EffectiveChargeDensity{T, N, S}, i::Int) where {T, N, S} = getindex(ρ.data, i)
@inline getindex(ρ::EffectiveChargeDensity{T, N, S}, s::Symbol) where {T, N, S} = getindex(ρ.grid, s)



function Base.NamedTuple(ρ::EffectiveChargeDensity{T}) where {T <: SSDFloat}
    # Legacy: the effective charge density is stored in units of V instead of the correct unit V/m²
    return (
        grid = NamedTuple(ρ.grid),
        values = ρ.data * internal_voltage_unit,
    )
end
Base.convert(T::Type{NamedTuple}, x::EffectiveChargeDensity) = T(x)
    
function EffectiveChargeDensity(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(ustrip(nt.values[1]))
    S = get_coordinate_system(grid)
    N = get_number_of_dimensions(grid)
    # Legacy: the effective charge density is parsed in units of V instead of the correct unit V/m²
    EffectiveChargeDensity{T, N, S, typeof(grid.axes)}( ustrip.(uconvert.(internal_voltage_unit, nt.values)), grid)
end
Base.convert(T::Type{EffectiveChargeDensity}, x::NamedTuple) = T(x)
