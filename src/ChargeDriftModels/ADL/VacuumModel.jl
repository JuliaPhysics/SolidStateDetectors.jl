mutable struct VacuumModel{T <: SSDFloat} <: AbstractTemperatureModel{T} end

function VacuumModel(config_file::Dict; T::Type{<:AbstractFloat} = Float32)::VacuumModel
    return VacuumModel{T}(config_file::Dict)
end

function VacuumModel{T}(config_file::Dict; temperature::Union{Missing, T} = missing)::VacuumModel where {T <: SSDFloat}
    m = VacuumModel{T}()
end

function scale_to_given_temperature(m::VacuumModel{T})::NTuple{4,T} where {T <: SSDFloat}
    return T(1), T(1), T(1), T(1)
end
