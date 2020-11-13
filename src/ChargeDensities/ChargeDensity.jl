abstract type AbstractChargeDensity{T <: SSDFloat} end

@inline function ChargeDensity(T::DataType, dict::Union{Dict{String, Any}, Dict{Any, Any}}, inputunit_dict::Dict)
    return ChargeDensity(T, Val{Symbol(dict["name"])}(), dict, inputunit_dict)
end

include("ConstantChargeDensity.jl")
include("LinearChargeDensity.jl")
include("CylindricalChargeDensity.jl")

