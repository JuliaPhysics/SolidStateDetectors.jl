struct DielectricDistribution{T, N, S} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S}
end

@inline size(ϵ::DielectricDistribution{T, N, S}) where {T, N, S} = size(ϵ.data)
@inline length(ϵ::DielectricDistribution{T, N, S}) where {T, N, S} = length(ϵ.data)
@inline getindex(ϵ::DielectricDistribution{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(ϵ.data, I...)
@inline getindex(ϵ::DielectricDistribution{T, N, S}, i::Int) where {T, N, S} = getindex(ϵ.data, i)
@inline getindex(ϵ::DielectricDistribution{T, N, S}, s::Symbol) where {T, N, S} = getindex(ϵ.grid, s)


function DielectricDistribution(fss::PotentialSimulationSetup{T, N, S})::DielectricDistribution{T, N, S} where {T, N, S}
    return DielectricDistribution{T, N, S}( fss.ϵ_r, fss.grid )
end


function NamedTuple(ϵ::DielectricDistribution{T}) where {T <: SSDFloat}
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
    DielectricDistribution{T, N, S}( nt.values, grid )
end
Base.convert(T::Type{DielectricDistribution}, x::NamedTuple) = T(x)
