"""
    ElectricFieldChargeDriftModel <: AbstractChargeDriftModel
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