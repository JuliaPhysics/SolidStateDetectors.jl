"""
    struct ConstantImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

Impurity density model that assumes a constant impurity density everywhere.

## Fields
* `ρ::T`: the constant value of the impurity density.

## Definition in Configuration File

A `ConstantImpurityDensity` is defined in the configuration file through the field `impurity_density` 
(of a `semiconductor`) with `name: constant` and a `value` for `ρ`.

An example definition of a constant impurity density looks like this:
```yaml 
impurity_density:
  name: constant
  value: 1.0e10 # 1/m³
```
"""
struct ConstantImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T} 
    ρ::T
end

function get_impurity_density(cdm::ConstantImpurityDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    return cdm.ρ
end

function ImpurityDensity(T::DataType, t::Val{:constant}, dict::AbstractDict, input_units::NamedTuple)
    ρ::T = haskey(dict, "value") ? _parse_value(T, dict["value"], input_units.length^(-3)) : T(0)
    ConstantImpurityDensity{T}( ρ )
end

(*)(scale::Real, lcdm::ConstantImpurityDensity{T}) where {T} = ConstantImpurityDensity{T}(T(scale * lcdm.ρ))

(+)(offset::Real, lcdm::ConstantImpurityDensity{T}) where {T} = ConstantImpurityDensity{T}(T(offset + lcdm.ρ))