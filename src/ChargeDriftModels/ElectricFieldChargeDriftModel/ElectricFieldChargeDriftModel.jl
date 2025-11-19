"""
    struct ElectricFieldChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
        
Charge drift model in which the electrons and holes drift along the electric field with a mobility of ± 1m²/Vs.

This model is the default when no charge drift model is defined in the configuration file.
"""
struct ElectricFieldChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T} end

function ElectricFieldChargeDriftModel(
        config::AbstractDict = Dict(), 
        input_units::Union{Missing, NamedTuple} = missing; 
        T::Type=Float32,
        temperature::Union{Missing, Real, Unitful.Temperature} = missing
    )::ElectricFieldChargeDriftModel{T}
   ElectricFieldChargeDriftModel{T}()
end

ElectricFieldChargeDriftModel{T}(args...; kwargs...) where {T <: SSDFloat} = ElectricFieldChargeDriftModel(args...; T=T, kwargs...)

getVe(fv::SVector{3, T}, cdm::ElectricFieldChargeDriftModel, current_pos::CartesianPoint{T} = zero(CartesianPoint{T})) where {T <: SSDFloat} = -fv
getVh(fv::SVector{3, T}, cdm::ElectricFieldChargeDriftModel, current_pos::CartesianPoint{T} = zero(CartesianPoint{T})) where {T <: SSDFloat} = fv