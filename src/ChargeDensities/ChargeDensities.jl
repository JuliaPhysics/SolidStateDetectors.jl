"""
    abstract type AbstractChargeDensity{T <: SSDFloat} end
    
Struct defining the charge density inside a [`Passive`](@ref).

For each charge density, there should be a method [`get_charge_density`](@ref) 
which returns the charge density in SI units (C/m³) at a given point `pt`.

## Examples
* [`ConstantChargeDensity`](@ref)
* [`LinearChargeDensity`](@ref)
* [`CylindricalChargeDensity`](@ref)
"""
abstract type AbstractChargeDensity{T <: SSDFloat} end

@inline function ChargeDensity(T::DataType, dict::AbstractDict, input_units::NamedTuple)
    update_charge_density_config!(dict, input_units)
    return ChargeDensity(T, Val{Symbol(dict["name"])}(), dict, input_units)
end

"""
    get_charge_density(cd::AbstractChargeDensity, pt::AbstractCoordinatePoint)
    
Returns the charge density at a given point, `pt`, based on the charge density model `cd`.

## Arguments
* `cd::AbstractChargeDensity`: The [`AbstractChargeDensity`](@ref) defining the charge density inside a [`Passive`](@ref).
* `pt::AbstractCoordinatePoint`: The point at which `cd` is to be evaluated.

!!! note 
    The value returned by `get_charge_density` is in units of C/m³ (SI units).
"""
function get_charge_density end

include("ConstantChargeDensity.jl")
include("LinearChargeDensity.jl")
include("CylindricalChargeDensity.jl")

