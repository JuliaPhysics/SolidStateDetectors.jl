"""
    struct CylindricalImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

Impurity density model which assumes a linear gradient in impurity density in each spatial dimension 
of a cylindrical coordinate system.
 
## Fields
* `offset::T`: impurity density value at the origin of the coordinate system.
* `gradients::NTuple{3,T}`: linear slopes in `r` and `z` direction.

## Definition in Configuration File

A `CylindricalImpurityDensity` is defined in the configuration file through the field `impurity_density` 
(of a `passive` or `surrounding`) with `name: cylindrical` and optional fields `r` and `z`
that can each contain `init` for initial values at 0 and `gradient` for gradients in that dimension.

An example definition of a cylindrical impurity density looks like this:
```yaml 
impurity_density:
  name: cylindrical
  r:  # impurity profile with linear gradient in r
    init: 1.0e10     # 1/m³
    gradient: 1.0e11 # 1/m⁴
```
"""
struct CylindricalImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}
    offset::T
    gradients::NTuple{3, T}
end

function ImpurityDensity(T::DataType, t::Val{:cylindrical}, dict::AbstractDict, input_units::NamedTuple)
    offset, gradients = zero(T), zeros(T,3)
    density_unit = input_units.length^(-3)
    density_gradient_unit = input_units.length^(-4)
    if prod(map(k -> k in ["r","z"], collect(keys(dict)))) @warn "Only r and z are supported in the cylindrical impurity density model.\nChange the impurity density model in the config file or remove all other entries." end
    if haskey(dict, "r")
        if haskey(dict["r"], "init")     offset += _parse_value(T, dict["r"]["init"], density_unit) end
        if haskey(dict["r"], "gradient") gradients[1] = _parse_value(T, dict["r"]["gradient"], density_gradient_unit) end
    end
    if haskey(dict, "z")
        if haskey(dict["z"], "init")     offset += _parse_value(T, dict["z"]["init"], density_unit) end
        if haskey(dict["z"], "gradient") gradients[3] = _parse_value(T, dict["z"]["gradient"], density_gradient_unit) end
    end
    CylindricalImpurityDensity{T}( T(offset), NTuple{3, T}(gradients) )
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