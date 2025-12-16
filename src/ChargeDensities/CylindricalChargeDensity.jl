"""
    struct CylindricalChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}

Charge density model which assumes a linear gradient in charge density in each spatial dimension 
of a cylindrical coordinate system.
 
## Fields
* `offset::T`: charge density values at the origin of the coordinate system.
* `gradients::NTuple{3,T}`: linear slopes in `r` and `z` direction.

## Definition in Configuration File

A `CylindricalChargeDensity` is defined in the configuration file through the field `charge_density` 
(of a `passive` or `surrounding`) with `name: cylindrical` and optional fields `offset` for the charge
density value at the origin of the coordinate system, and `gradient` for gradients in `r` and/or `z`.

An example definition of a cylindrical charge density looks like this:
```yaml 
# charge density with linear gradient in r
charge_density:
  name: cylindrical
  offset: 1.0e-10     # C/m³
  gradient: 
    r: 1.0e-11 # C/m⁴
```
"""
struct CylindricalChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}
    offset::T
    gradients::NTuple{3, T}
end

function ChargeDensity(T::DataType, t::Val{:cylindrical}, dict::AbstractDict, input_units::NamedTuple)
    density_unit = internal_charge_unit * input_units.length^(-3)
    density_gradient_unit = internal_charge_unit * input_units.length^(-4)
    offset::T = _parse_value(T, get(dict, "offset", get(dict, "init", 0)), density_unit)
    gradients::NTuple{3, T} = let g = get(dict, "gradient", get(dict, "gradients", Dict()))
        haskey(g, "phi") && @warn "Ignoring gradient for phi in cylindrical charge density"
        _parse_value.(T, (get(g, "r", 0), 0, get(g, "z", 0)), density_gradient_unit)
    end
    CylindricalChargeDensity{T}( offset, gradients )
end

function get_charge_density(lcdm::CylindricalChargeDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    pt::CylindricalPoint{T} = CylindricalPoint(pt)
    ρ::T = lcdm.offset
    for i in eachindex(lcdm.gradients)
        ρ += pt[i] * lcdm.gradients[i]
    end
    return ρ
end
