"""
    struct CylindricalChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}

Charge density model which assumes a linear gradient in charge density in each spatial dimension 
of a cylindrical coordinate system.
 
## Fields
* `offsets::NTuple{3,T}`: charge density values at the origin of each dimension.
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
    offsets::NTuple{3, T}
    gradients::NTuple{3, T}
end

function ChargeDensity(T::DataType, t::Val{:cylindrical}, dict::AbstractDict, input_units::NamedTuple)
    offsets, gradients = zeros(T,3), zeros(T,3)
    density_unit = internal_charge_density * input_units.length^(-3)
    density_gradient_unit = internal_charge_density * input_units.length^(-4)
    if prod(map(k -> k in ["r","z"], collect(keys(dict)))) @warn "Only r and z are supported in the cylindrical charge density model.\nChange the charge density model in the config file or remove all other entries." end
    if haskey(dict, "r")     
        if haskey(dict["r"], "init")     offsets[1]   = _parse_value(T, dict["r"]["init"], density_unit) end
        if haskey(dict["r"], "gradient") gradients[1] = _parse_value(T, dict["r"]["gradient"], density_gradient_unit) end
    end
    if haskey(dict, "z")     
        if haskey(dict["z"], "init")     offsets[3]   = _parse_value(T, dict["z"]["init"], density_unit) end
        if haskey(dict["z"], "gradient") gradients[3] = _parse_value(T, dict["z"]["gradient"], density_gradient_unit) end
    end
    CylindricalChargeDensity{T}( NTuple{3, T}(offsets), NTuple{3, T}(gradients) )
end

function get_charge_density(lcdm::CylindricalChargeDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    pt::CylindricalPoint{T} = CylindricalPoint(pt)
    ρ::T = 0
    for i in eachindex(lcdm.offsets)
        ρ += (lcdm.offsets[i] + pt[i] * lcdm.gradients[i]) #* T(1e16) # * T(1e10) * T(1e6) -> 1/cm^3 -> 1/m^3
    end
    return ρ
end
