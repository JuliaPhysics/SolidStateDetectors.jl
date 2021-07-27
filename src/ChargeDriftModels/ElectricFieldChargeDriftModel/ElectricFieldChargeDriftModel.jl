"""
    struct ElectricFieldChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
        
Charge drift model in which the electrons and holes drift along the electric field with a mobility of ± 1m²/Vs.

This model is the default when no charge drift model is defined in the configuration file.
"""
struct ElectricFieldChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T} end

ElectricFieldChargeDriftModel(T::Type{<:SSDFloat}) = ElectricFieldChargeDriftModel{T}()

function get_electron_drift_field(ef::Array{SVector{3,T},3}, ::ElectricFieldChargeDriftModel)::Array{SVector{3,T},3} where {T <: SSDFloat}
    return -ef
end
function get_hole_drift_field(ef::Array{SVector{3,T},3}, ::ElectricFieldChargeDriftModel)::Array{SVector{3,T},3} where {T <: SSDFloat}
    return ef
end

function getVe(fv::SVector{3, T}, cdm::ElectricFieldChargeDriftModel)::SVector{3, T} where {T <: SSDFloat}
    return -fv
end
function getVh(fv::SVector{3, T}, cdm::ElectricFieldChargeDriftModel)::SVector{3, T} where {T <: SSDFloat}
    return fv
end