"""
    struct LinearImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

Impurity density model which assumes a linear gradient in impurity density in each dimension of a Cartesian coordinate system.
 
## Fields
* `offsets::NTuple{3,T}`: impurity density values at the origin of each dimension.
* `gradients::NTuple{3,T}`: linear slopes in `x`, `y` and `z` direction.

## Definition in Configuration File

A `LinearImpurityDensity` is defined in the configuration file through the field `impurity_density` 
(of a `passive` or `surrounding`) with `name: linear` and optional fields `x`, `y` and `z`
that can each contain `init` for initial values at 0 and `gradient` for gradients in that dimension.

An example definition of a linear impurity density looks like this:
```yaml 
impurity_density:
  name: linear
  x:  # impurity profile with linear gradient in x
    init: 1.0e10     # 1/m³
    gradient: 1.0e11 # 1/m⁴
```
"""
struct LinearImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}
    offsets::NTuple{3, T}
    gradients::NTuple{3, T}
end

function ImpurityDensity(T::DataType, t::Val{:linear}, dict::AbstractDict, input_units::NamedTuple)
    offsets, gradients = zeros(T,3), zeros(T,3)
    density_unit = input_units.length^(-3)
    density_gradient_unit = input_units.length^(-4)
    if prod(map(k -> k in ["x","y","z"], collect(keys(dict)))) @warn "Only x, y and z are supported in the linear impurity density model.\nChange the impurity density model in the config file or remove all other entries." end
    if haskey(dict, "x")     
        if haskey(dict["x"], "init")     offsets[1]   = _parse_value(T, dict["x"]["init"], density_unit) end
        if haskey(dict["x"], "gradient") gradients[1] = _parse_value(T, dict["x"]["gradient"], density_gradient_unit) end
    end
    if haskey(dict, "y")     
        if haskey(dict["y"], "init")     offsets[2]   = _parse_value(T, dict["y"]["init"], density_unit) end
        if haskey(dict["y"], "gradient") gradients[2] = _parse_value(T, dict["y"]["gradient"], density_gradient_unit) end
    end
    if haskey(dict, "z")     
        if haskey(dict["z"], "init")     offsets[3]   = _parse_value(T, dict["z"]["init"], density_unit) end
        if haskey(dict["z"], "gradient") gradients[3] = _parse_value(T, dict["z"]["gradient"], density_gradient_unit) end
    end
    LinearImpurityDensity{T}( NTuple{3, T}(offsets), NTuple{3, T}(gradients) )
end

function get_impurity_density(lcdm::LinearImpurityDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    pt::CartesianPoint{T} = CartesianPoint(pt)
    ρ::T = 0
    for i in eachindex(lcdm.offsets)
        ρ += (lcdm.offsets[i] + pt[i] * lcdm.gradients[i]) #* T(1e16) # * T(1e10) * T(1e6) -> 1/cm^3 -> 1/m^3
    end
    return ρ
end

(*)(scale::Real, lcdm::LinearImpurityDensity{T}) where {T} = LinearImpurityDensity{T}(T.(scale .* lcdm.offsets), T.(scale .* lcdm.gradients))

(+)(offset::Union{<:Real, <:Quantity{<:Real, Unitful.𝐋^(-3)}}, lcdm::LinearImpurityDensity{T}) where {T} = LinearImpurityDensity{T}(T.(to_internal_units(offset) .+ lcdm.offsets), lcdm.gradients)
