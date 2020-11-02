struct ChargeDensity{T, N, S} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S}
end

@inline size(ρ::ChargeDensity{T, N, S}) where {T, N, S} = size(ρ.data)
@inline length(ρ::ChargeDensity{T, N, S}) where {T, N, S} = length(ρ.data)
@inline getindex(ρ::ChargeDensity{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(ρ.data, I...)
@inline getindex(ρ::ChargeDensity{T, N, S}, i::Int) where {T, N, S} = getindex(ρ.data, i)
@inline getindex(ρ::ChargeDensity{T, N, S}, s::Symbol) where {T, N, S} = getindex(ρ.grid, s)


function ChargeDensity(fss::PotentialSimulationSetup{T, N, S})::ChargeDensity{T, N, S} where {T, N, S}
    return ChargeDensity{T, N, S}( fss.ρ, fss.grid )
end


function NamedTuple(ρ::ChargeDensity{T}) where {T <: SSDFloat}
    return (
        grid = NamedTuple(ρ.grid),
        values = ρ.data * internal_voltage_unit,
    )
end
Base.convert(T::Type{NamedTuple}, x::ChargeDensity) = T(x)
    
function ChargeDensity(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(ustrip(nt.values[1]))
    S = get_coordinate_system(grid)
    N = get_number_of_dimensions(grid)
    ChargeDensity{T, N, S}( ustrip.(uconvert.(internal_voltage_unit, nt.values)), grid)
end
Base.convert(T::Type{ChargeDensity}, x::NamedTuple) = T(x)
