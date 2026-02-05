"""
    struct LinearChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}

Charge density model which assumes a linear gradient in charge density in each dimension of a Cartesian coordinate system.
 
## Fields
* `offset::T`: charge density values at the origin of the coordinate system.
* `gradients::NTuple{3,T}`: linear slopes in `x`, `y` and `z` direction.

## Definition in Configuration File

A `LinearChargeDensity` is defined in the configuration file through the field `charge_density` 
(of a `passive` or `surrounding`) with `name: linear` and optional fields `offset` for the charge
density value at the origin of the coordinate system, and `gradient` for gradients in `x`, `y`,` and/or `z`.

An example definition of a linear charge density looks like this:
```yaml 
# charge density with linear gradient in x
charge_density:
  name: linear
  offset: 1.0e-10 # C/m³
  gradient: 
    x: 1.0e-11  # C/m⁴
```
"""
struct LinearChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}
    offset::T
    gradients::NTuple{3, T}
end

function ChargeDensity(T::DataType, t::Val{:linear}, dict::AbstractDict, input_units::NamedTuple)
    density_unit = internal_charge_unit * input_units.length^(-3)
    density_gradient_unit = internal_charge_unit * input_units.length^(-4)
    offset::T = _parse_value(T, get(dict, "offset", get(dict, "init", 0)), density_unit)
    gradients::NTuple{3, T} = let g = get(dict, "gradient", get(dict, "gradients", Dict()))
        _parse_value.(T, (get(g, "x", 0), get(g, "y", 0), get(g, "z", 0)), density_gradient_unit)
    end
    LinearChargeDensity{T}( offset, gradients )
end

function get_charge_density(lcdm::LinearChargeDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    pt::CartesianPoint{T} = CartesianPoint(pt)
    ρ::T = lcdm.offset
    for i in eachindex(lcdm.gradients)
        ρ += pt[i] * lcdm.gradients[i]
    end
    return ρ
end
