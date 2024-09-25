"""
    struct ElectricFieldChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
        
Charge drift model in which the electrons and holes drift along the electric field with a mobility of ± 1m²/Vs.

This model is the default when no charge drift model is defined in the configuration file.
"""
struct ElectricFieldChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T} 
    ElectricFieldChargeDriftModel{T}(::AbstractDict = Dict()) where {T <: SSDFloat} = new{T}()
end

ElectricFieldChargeDriftModel(T::Type{<:SSDFloat}) = ElectricFieldChargeDriftModel{T}()
getVe(fv::SVector{3, T}, cdm::ElectricFieldChargeDriftModel) where {T <: SSDFloat} = -fv
getVh(fv::SVector{3, T}, cdm::ElectricFieldChargeDriftModel) where {T <: SSDFloat} = fv