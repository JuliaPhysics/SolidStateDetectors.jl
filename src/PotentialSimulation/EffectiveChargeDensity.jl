"""
    struct EffectiveChargeDensity{T, N, S, AT} <: AbstractArray{T, N}
- `T`: Element type of `data`.
- `N`: Dimension of the `grid` and `data` array.  
- `S`: Coordinate system (`Cartesian` or `Cylindrical`).
- `AT`: Axes type.  
        
# Fields
- `data::Array{T, N}`
- `grid::Grid{T, N, S, AT}`

The `data` array contains the values of the effective charge density at the discrete points defined by the 
axes ticks of the `grid`.

The effective charge density is the charge density ([Q/m^3]) times the volume of the voxel of the respective
grid point ([m^3]). Thus, the units of the effective charge density in `data` is coulombs ([C]).
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


function EffectiveChargeDensity(fss::PotentialSimulationSetup{T, N, S})::EffectiveChargeDensity{T, N, S} where {T, N, S}
    return EffectiveChargeDensity{T, N, S}( fss.ρ, fss.grid )
end


function NamedTuple(ρ::EffectiveChargeDensity{T}) where {T <: SSDFloat}
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
    EffectiveChargeDensity{T, N, S}( ustrip.(uconvert.(internal_voltage_unit, nt.values)), grid)
end
Base.convert(T::Type{EffectiveChargeDensity}, x::NamedTuple) = T(x)
