mutable struct VacuumModel{T <: SSDFloat} <: AbstractTemperatureModel{T} end

function VacuumModel(config_file::AbstractDict; T::Type{<:AbstractFloat} = Float32)::VacuumModel
    return VacuumModel{T}(config_file::AbstractDict)
end

function VacuumModel{T}(config_file::AbstractDict)::VacuumModel where {T <: SSDFloat}
    m = VacuumModel{T}()
end

function scale_to_temperature(cdm::ADL2016ChargeDriftModel{T,M,N,VacuumModel{T}}, Temp::RealQuantity) where {T,M,N}
    cdm
end