mutable struct VacuumModel{T<:AbstractFloat} <: TemperatureModels{T} end

function VacuumModel(config_file::Dict; T::Type{<:AbstractFloat} = Float32)::VacuumModel
    return VacuumModel{T}(config_file::Dict)
end

function VacuumModel{T}(config_file::Dict; temperature::Union{Missing, T} = missing)::VacuumModel where T<:AbstractFloat
    m = VacuumModel{T}()
end

function scale_to_given_temperature(m::VacuumModel{T})::NTuple{4,T} where T<:AbstractFloat
    return 1.0,1.0,1.0,1.0
end
