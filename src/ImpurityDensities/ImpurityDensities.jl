"""
    abstract type AbstractImpurityDensity{T <: SSDFloat} end
    
Struct defining an impurity density inside a [`Semiconductor`](@ref).

For each impurity density, there should be a method
[`get_impurity_density`](@ref) which returns the impurity density in SI units (1/m³) at a given point.

## Examples
* [`ConstantImpurityDensity`](@ref)
* [`LinearImpurityDensity`](@ref)
* [`CylindricalImpurityDensity`](@ref)

Boule impurity densities are meant to be used when the impurity density is defined in the boule coordinates, where the z-axis is aligned with the boule growth direction. 

Different models are provided. In each the field `det_z0` is the z-coordinate of the detector origin in boule coordinates. The z-direction of the detector is opposite to the z-direction of the boule coordinates.

In this matter the detector impurities are automatically determined from those of the boule, depending on `det_z0`.

!!! note
    The sign of the impurity density is important. It is taken into account in the conversion to a charge density
    and, thus, defines where the semiconductor is n-type or p-type. 
"""
abstract type AbstractImpurityDensity{T <: SSDFloat} end

@inline function ImpurityDensity(T::DataType, dict::AbstractDict, input_units::NamedTuple)
    return ImpurityDensity(T, Val{Symbol(dict["name"])}(), dict, input_units)
end

"""
    get_impurity_density(id::AbstractImpurityDensity, pt::AbstractCoordinatePoint)
    
Returns the impurity density at a given point, `pt`, based on the impurity density model `id`.

## Arguments
* `id::AbstractImpurityDensity`: The [`AbstractImpurityDensity`](@ref) defining the impurity density inside a [`Semiconductor`](@ref).
* `pt::AbstractCoordinatePoint`: The point at which `id` is to be evaluated.

!!! note 
    The value returned by `get_impurity_density` is in units of 1/m³ (SI Units).
"""
function get_impurity_density end


include("ConstantImpurityDensity.jl")
include("LinearImpurityDensity.jl")
include("CylindricalImpurityDensity.jl")

include("BouleImpurityDensities/BouleImpurityDensities.jl")