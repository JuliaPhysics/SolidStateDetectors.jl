"""
    # abstract type AbstractChargeDensity{T <: SSDFloat} end

Charge densities <: `AbstractChargeDensity` should be defined to return (via `get_charge_density`-method)
a charge density in SI units, thus, in Q/m^3.
"""
abstract type AbstractChargeDensity{T <: SSDFloat} end

@inline function ChargeDensity(T::DataType, dict::Union{Dict{String, Any}, Dict{Any, Any}}, input_units::NamedTuple)
    return ChargeDensity(T, Val{Symbol(dict["name"])}(), dict, input_units)
end

include("ConstantChargeDensity.jl")
include("LinearChargeDensity.jl")
include("CylindricalChargeDensity.jl")

