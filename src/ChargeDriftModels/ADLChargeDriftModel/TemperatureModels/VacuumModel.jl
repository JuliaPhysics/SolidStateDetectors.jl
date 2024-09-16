mutable struct VacuumModel{T <: SSDFloat} <: AbstractTemperatureModel{T} end

function VacuumModel(config_file::AbstractDict; T::Type{<:AbstractFloat} = Float32)::VacuumModel
    return VacuumModel{T}(config_file::AbstractDict)
end

function VacuumModel{T}(config_file::AbstractDict; temperature::Union{Missing, T} = missing)::VacuumModel where {T <: SSDFloat}
    m = VacuumModel{T}()
end

function scale_to_given_temperature(m::VacuumModel{T})::NTuple{4,T} where {T <: SSDFloat}
    return T(1), T(1), T(1), T(1)
end

print(io::IO, tm::VacuumModel{T}) where {T <: SSDFloat} = print(io, "No temperature model defined")
println(io::IO, tm::VacuumModel) = print(io, tm)
