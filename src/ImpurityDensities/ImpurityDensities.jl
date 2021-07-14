"""
    # abstract type AbstractImpurityDensity{T <: SSDFloat} end

Impurity densities <: `AbstractImpurityDensity` should be defined to return (via `get_impurity_density`-method)
a density of impurities in SI units, thus, in 1/m^3. For semiconductors, this will be converted internally
into a charge distribution by multiplication with the elementary charge. 
The sign of the impurity density is important. It is taken into account in the conversion to a charge density
and, thus, defines where the semiconductor is n-type or p-type. 
"""
abstract type AbstractImpurityDensity{T <: SSDFloat} end

@inline function ImpurityDensity(T::DataType, dict::Union{Dict{String, Any}, Dict{Any, Any}}, input_units::NamedTuple)
    return ImpurityDensity(T, Val{Symbol(dict["name"])}(), dict, input_units)
end

include("ConstantImpurityDensity.jl")
include("LinearImpurityDensity.jl")
include("CylindricalImpurityDensity.jl")

