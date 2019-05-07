abstract type AbstractChargeDriftModel{T <: SSDFloat} end
abstract type AbstractTemperatureModel{T <: SSDFloat} end


include("Vacuum/Vacuum.jl")
include("ADL/ADL.jl")