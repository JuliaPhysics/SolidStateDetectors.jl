"""
    struct CylindricalImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

Impurity density model which assumes a linear gradient in impurity density in each spatial dimension 
of a cylindrical coordinate system.
 
## Fields
* `offset::T`: impurity density value at the origin of the coordinate system.
* `gradients::NTuple{3,T}`: linear slopes in `r` and `z` direction.

## Definition in Configuration File

A `CylindricalImpurityDensity` is defined in the configuration file through the field `impurity_density` 
(of a `passive` or `surrounding`) with `name: cylindrical`  and optional fields `offset` for the impurity
density value at the origin of the coordinate system, and `gradient` for gradients in `r` and/or `z`.

An example definition of a cylindrical impurity density looks like this:
```yaml
# impurity profile with linear gradient in r 
impurity_density:
  name: cylindrical
  offset: 1.0e10 # 1/m³
  gradient: 
    r: 1.0e11  # 1/m⁴
```
"""
struct CylindricalImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}
    offset::T
    gradients::NTuple{3, T}
end

function ImpurityDensity(T::DataType, t::Val{:cylindrical}, dict::AbstractDict, input_units::NamedTuple)
    density_unit = input_units.length^(-3)
    density_gradient_unit = input_units.length^(-4)
    offset::T = _parse_value(T, get(dict, "offset", get(dict, "init", 0)), density_unit)
    gradients::NTuple{3, T} = let g = get(dict, "gradient", get(dict, "gradients", Dict()))
        haskey(g, "phi") && @warn "Ignoring gradient for phi in cylindrical impurity density"
        _parse_value.(T, (get(g, "r", 0), 0, get(g, "z", 0)), density_gradient_unit)
    end
    CylindricalImpurityDensity{T}( offset, gradients )
end

function get_impurity_density(lcdm::CylindricalImpurityDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    pt::CylindricalPoint{T} = CylindricalPoint(pt)
    ρ::T = lcdm.offset
    for i in eachindex(lcdm.gradients)
        ρ += pt[i] * lcdm.gradients[i]
    end
    return ρ
end

(*)(scale::Real, lcdm::CylindricalImpurityDensity{T}) where {T} = CylindricalImpurityDensity{T}(T(scale * lcdm.offset), T.(scale .* lcdm.gradients))

(+)(offset::Union{<:Real, <:Quantity{<:Real, Unitful.𝐋^(-3)}}, lcdm::CylindricalImpurityDensity{T}) where {T} = CylindricalImpurityDensity{T}(T(to_internal_units(offset) + lcdm.offset), lcdm.gradients)