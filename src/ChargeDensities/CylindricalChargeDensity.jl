"""
    struct CylindricalChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}

Charge density model which assumes a linear gradient in charge density in each spatial dimension 
of a cylindrical coordinate system.
 
## Fields
* `offset::T`: charge density values at the origin of the coordinate system.
* `gradients::NTuple{3,T}`: linear slopes in `r` and `z` direction.

## Definition in Configuration File

A `CylindricalChargeDensity` is defined in the configuration file through the field `charge_density` 
(of a `passive` or `surrounding`) with `name: cylindrical` and optional fields `r` and `z`
that can each contain `init` for initial values at 0 and `gradient` for gradients in that dimension.

An example definition of a cylindrical charge density looks like this:
```yaml 
charge_density:
  name: cylindrical
  r:  # charge density with linear gradient in r
    init: 1.0e-10     # C/m³
    gradient: 1.0e-11 # C/m⁴
```
"""
struct CylindricalChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}
    offset::T
    gradients::NTuple{3, T}
end

function ChargeDensity(T::DataType, t::Val{:cylindrical}, dict::AbstractDict, input_units::NamedTuple)
    offset, gradients = zero(T), zeros(T,3)
    density_unit = internal_charge_unit * input_units.length^(-3)
    density_gradient_unit = internal_charge_unit * input_units.length^(-4)
    if prod(map(k -> k in ["r","z"], collect(keys(dict)))) @warn "Only r and z are supported in the cylindrical charge density model.\nChange the charge density model in the config file or remove all other entries." end
    if haskey(dict, "r")
        if haskey(dict["r"], "init")     offset += _parse_value(T, dict["r"]["init"], density_unit) end
        if haskey(dict["r"], "gradient") gradients[1] = _parse_value(T, dict["r"]["gradient"], density_gradient_unit) end
    end
    if haskey(dict, "z")     
        if haskey(dict["z"], "init")     offset +=  _parse_value(T, dict["z"]["init"], density_unit) end
        if haskey(dict["z"], "gradient") gradients[3] = _parse_value(T, dict["z"]["gradient"], density_gradient_unit) end
    end
    CylindricalChargeDensity{T}( T(offset), NTuple{3, T}(gradients) )
end

function get_charge_density(lcdm::CylindricalChargeDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    pt::CylindricalPoint{T} = CylindricalPoint(pt)
    ρ::T = lcdm.offset
    for i in eachindex(lcdm.gradients)
        ρ += pt[i] * lcdm.gradients[i]
    end
    return ρ
end
