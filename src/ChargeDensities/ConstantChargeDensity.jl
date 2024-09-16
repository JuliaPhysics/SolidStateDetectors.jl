"""
    struct ConstantChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}

Charge density model that assumes a constant charge density everywhere.

## Fields
* `ρ::T`: the constant value of the charge density.

## Definition in Configuration File

A `ConstantChargeDensity` is defined in the configuration file through the field `charge_density` 
(of a `passive` or `surrounding`) with `name: constant` and a `value` for `ρ`.

An example definition of a constant charge density looks like this:
```yaml 
charge_density:
  name: constant
  value: 1.0e-10 # C/m³
```
"""
struct ConstantChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T} 
    ρ::T
end

function get_charge_density(cdm::ConstantChargeDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    return cdm.ρ
end


function ChargeDensity(T::DataType, t::Val{:constant}, dict::AbstractDict, input_units::NamedTuple)
    ρ::T = haskey(dict, "value") ? _parse_value(T, dict["value"], input_units.length^(-3)) : T(0)
    ConstantChargeDensity{T}( ρ )
end
