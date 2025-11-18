mutable struct VacuumTemperatureModel{T <: SSDFloat} <: AbstractTemperatureModel{T} end

function VacuumTemperatureModel(config::AbstractDict; T::Type{<:AbstractFloat} = Float32)::VacuumTemperatureModel
    VacuumTemperatureModel{T}(config)
end

function VacuumTemperatureModel{T}(::AbstractDict)::VacuumTemperatureModel where {T <: SSDFloat}
    VacuumTemperatureModel{T}()
end

function TemperatureModel(T::DataType, ::Val{:Vacuum}, ::AbstractDict)
    VacuumTemperatureModel{T}()
end

function scale_to_temperature(cdm::CDM, Temp::Union{<:Real,Unitful.Temperature}) where {T,M,N,CDM<:Union{ADL2016ChargeDriftModel{T,M,N,VacuumTemperatureModel{T}},ADLChargeDriftModel{T,M,N,VacuumTemperatureModel{T}}}}
    return cdm
end