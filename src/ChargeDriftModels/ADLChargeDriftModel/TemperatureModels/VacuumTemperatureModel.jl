mutable struct VacuumTemperatureModel{T <: SSDFloat} <: AbstractTemperatureModel{T} end

function VacuumTemperatureModel(::AbstractDict = Dict(), ::Union{Missing, NamedTuple} = missing; T::Type{<:AbstractFloat} = Float32)
    VacuumTemperatureModel{T}()
end

VacuumTemperatureModel{T}(args...; kwargs...) where {T <: SSDFloat} = VacuumTemperatureModel(args...; T=T, kwargs...)

function TemperatureModel(T::DataType, ::Val{:Vacuum}, ::AbstractDict, ::Union{Missing, NamedTuple})
    VacuumTemperatureModel{T}()
end

function scale_to_temperature(cdm::CDM, ::Union{<:Real,Unitful.Temperature}) where {T,M,N,CDM<:Union{ADL2016ChargeDriftModel{T,M,N,VacuumTemperatureModel{T}},ADLChargeDriftModel{T,M,N,VacuumTemperatureModel{T}}}}
    cdm
end