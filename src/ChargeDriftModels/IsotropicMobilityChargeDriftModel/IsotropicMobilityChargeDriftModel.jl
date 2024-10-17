"""
    struct IsotropicMobilityChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
        
Charge drift model in which the electrons and holes drift along the electric field with a constant mobility in m²/Vs.

"""
struct IsotropicMobilityChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
    μ_e::T
    μ_h::T
    IsotropicMobilityChargeDriftModel{T}(μ_e::RealQuantity, μ_h::RealQuantity) where {T <: SSDFloat} = new{T}(
        T(ustrip(u"m^2/V/s", μ_e)),
        T(ustrip(u"m^2/V/s", μ_h))
    )
end

@fastmath function getVe(fv::SVector{3, T}, cdm::IsotropicMobilityChargeDriftModel{T}) where {T <: SSDFloat}
    @inbounds begin 
        -cdm.μ_e*fv
    end
end

@fastmath function getVh(fv::SVector{3, T}, cdm::IsotropicMobilityChargeDriftModel{T}) where {T <: SSDFloat}
    @inbounds begin 
        cdm.μ_h*fv
    end
end