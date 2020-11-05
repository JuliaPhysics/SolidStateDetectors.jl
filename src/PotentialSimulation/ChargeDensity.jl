struct EffectiveChargeDensity{T, N, S} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S}
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
