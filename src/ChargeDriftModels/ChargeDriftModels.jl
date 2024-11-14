abstract type AbstractChargeDriftModel{T <: SSDFloat} end
abstract type AbstractTemperatureModel{T <: SSDFloat} end

include("ElectricFieldChargeDriftModel/ElectricFieldChargeDriftModel.jl")
include("ADLChargeDriftModel/ADLChargeDriftModel.jl")
include("IsotropicChargeDriftModel/IsotropicChargeDriftModel.jl")
