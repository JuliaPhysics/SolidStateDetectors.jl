"""
    struct LinearChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}

Charge density model which assumes a linear gradient in charge density in each dimension of a Cartesian coordinate system.
 
## Fields
* `offset::T`: charge density values at the origin of the coordinate system.
* `gradients::NTuple{3,T}`: linear slopes in `x`, `y` and `z` direction.

## Definition in Configuration File

A `LinearChargeDensity` is defined in the configuration file through the field `charge_density` 
(of a `passive` or `surrounding`) with `name: linear` and optional fields `x`, `y` and `z`
that can each contain `init` for initial values at 0 and `gradient` for gradients in that dimension.

An example definition of a linear charge density looks like this:
```yaml 
# charge density with linear gradient in x
charge_density:
  name: linear
  init: 1.0e-10 # C/m³
  gradient: 
    x: 1.0e-11  # C/m⁴
```
"""
struct LinearChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}
    offset::T
    gradients::NTuple{3, T}
end

function ChargeDensity(T::DataType, t::Val{:linear}, dict::AbstractDict, input_units::NamedTuple)
    offset, gradients = zero(T), zeros(T,3)
    density_unit = internal_charge_unit * input_units.length^(-3)
    density_gradient_unit = internal_charge_unit * input_units.length^(-4)
    if prod(map(k -> k in ["x","y","z"], collect(keys(dict)))) @warn "Only x, y and z are supported in the linear charge density model.\nChange the charge density model in the config file or remove all other entries." end
    if haskey(dict, "x")
        if haskey(dict["x"], "init")     offset += _parse_value(T, dict["x"]["init"], density_unit) end
        if haskey(dict["x"], "gradient") gradients[1] = _parse_value(T, dict["x"]["gradient"], density_gradient_unit) end
    end
    if haskey(dict, "y")
        if haskey(dict["y"], "init")     offset += _parse_value(T, dict["y"]["init"], density_unit) end
        if haskey(dict["y"], "gradient") gradients[2] = _parse_value(T, dict["y"]["gradient"], density_gradient_unit) end
    end
    if haskey(dict, "z")
        if haskey(dict["z"], "init")     offset += _parse_value(T, dict["z"]["init"], density_unit) end
        if haskey(dict["z"], "gradient") gradients[3] = _parse_value(T, dict["z"]["gradient"], density_gradient_unit) end
    end
    LinearChargeDensity{T}( T(offset), NTuple{3, T}(gradients) )
end

function get_charge_density(lcdm::LinearChargeDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    pt::CartesianPoint{T} = CartesianPoint(pt)
    ρ::T = lcdm.offset
    for i in eachindex(lcdm.gradients)
        ρ += pt[i] * lcdm.gradients[i]
    end
    return ρ
end
