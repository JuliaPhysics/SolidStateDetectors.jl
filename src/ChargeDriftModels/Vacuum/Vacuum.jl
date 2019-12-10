"""
    VacuumChargeDriftModel <: AbstractChargeDriftModel
"""
struct VacuumChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T} end

function get_electron_drift_field(ef::Array{SVector{3,T},3}, ::VacuumChargeDriftModel)::Array{SVector{3,T},3} where {T <: SSDFloat}
    return -ef
end
function get_hole_drift_field(ef::Array{SVector{3,T},3}, ::VacuumChargeDriftModel)::Array{SVector{3,T},3} where {T <: SSDFloat}
    return ef
end

function getVe(fv::SVector{3, T}, cdm::VacuumChargeDriftModel)::SVector{3, T} where {T <: SSDFloat}
    return -fv
end
function getVh(fv::SVector{3, T}, cdm::VacuumChargeDriftModel)::SVector{3, T} where {T <: SSDFloat}
    return fv
end