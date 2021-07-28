"""
    struct Box{T, CO} <: AbstractVolumePrimitive{T}

Volume primitive describing a three-dimensional box with its surfaces
being parallel to the `xy`, `xy` and `yz` plane.

## Parametric types
* `T`: Precision type.
* `CO`: Describes whether the surface belongs to the primitive. 
    It can be `ClosedPrimitive`, i.e. the surface points belong to the primitive,
    or `OpenPrimitive`, i.e. the surface points do not belong to the primitive.
    
## Fields
* `hX::T`: Half of the width in `x` dimension (in m).
* `hY::T`: Half of the width in `y` dimension (in m).
* `hZ::T`: Half of the width in `z` dimension (in m).
* `origin::CartesianPoint{T}`: The position of the center of the `Box`.
* `rotation::SMatrix{3,3,T,9}`: Matrix that describes a rotation of the `Box` around its `origin`.

## Definition in Configuration File

A `Box` is defined in the configuration file as part of the `geometry` field 
of an object through the field `box`.

Example definitions of a `Box` looks like this:
```yaml
box:
  widths: [2, 4, 6] # => hX = 1; hY = 2; hZ = 3;
  origin: [0, 0, 0] # [x, y, z] - Optional; Default: [0, 0, 0]
  rotate: # Optional; Default: no rotation
    Z: 0 
```

The halfwidths `hX`, `hY` and `hZ` can also be defined directly in the configuration file:
```yaml
box:
  halfwidths: [1, 2, 3] # => hX = 1; hY = 2; hZ = 3;
```
or
```yaml
box:
  hX: 1
  hY: 2
  hZ: 3
```
"""
@with_kw struct Box{T, CO} <: AbstractVolumePrimitive{T, CO}
    hX::T = 1
    hY::T = 1
    hZ::T = 1
    origin::CartesianPoint{T} = zero(CartesianPoint{T})
    rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})
end

Box{T, CO}( b::Box{T, CO}; COT = CO,
            origin::CartesianPoint{T} = b.origin,
            rotation::SMatrix{3,3,T,9} = b.rotation) where {T, CO<:Union{ClosedPrimitive, OpenPrimitive}} =
    Box{T, COT}(b.hX, b.hY, b.hZ, origin, rotation)

function _in(pt::CartesianPoint{T}, b::Box{<:Any, ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T}
    abs(pt.x) <= b.hX + csgtol && 
    abs(pt.y) <= b.hY + csgtol && 
    abs(pt.z) <= b.hZ + csgtol 
end

_in(pt::CartesianPoint{T}, b::Box{<:Any, OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} = 
    abs(pt.x) < b.hX - csgtol && abs(pt.y) < b.hY - csgtol && abs(pt.z) < b.hZ - csgtol
 

function Geometry(::Type{T}, ::Type{Box}, dict::AbstractDict, input_units::NamedTuple, transformations::Transformations{T}) where {T}
    length_unit = input_units.length
    origin = get_origin(T, dict, length_unit)
    rotation = get_rotation(T, dict, input_units.angle) 
    hX, hY, hZ = if haskey(dict, "widths")
        @assert length(dict["widths"]) == 3
        _parse_value(T, dict["widths"], length_unit)/2
    elseif haskey(dict, "halfwidths")
        @assert length(dict["halfwidths"]) == 3
        _parse_value(T, dict["halfwidths"], length_unit)
    elseif haskey(dict, "hX") && haskey(dict, "hZ") && haskey(dict, "hZ")
        _parse_value(T, dict["hX"], length_unit), 
        _parse_value(T, dict["hY"], length_unit), 
        _parse_value(T, dict["hZ"], length_unit)
    elseif haskey(dict, "x") && haskey(dict, "y") && haskey(dict, "z")
        @warn "Deprecation warning: Detected old primitive definition for `Box`. 
            Please update your configuration file to the new format 
            via `widths`, `halfwidths` or `hX`, `hY` and `hZ`
            in combination with possible fields `origin` and `rotate`.
            The old definition overwrites the optional field `origin`."        
        x = parse_interval_of_primitive(T, "x", dict, length_unit)
        y = parse_interval_of_primitive(T, "y", dict, length_unit)
        z = parse_interval_of_primitive(T, "z", dict, length_unit)
        μx = typeof(x) <: Real ? zero(T) : mean(x)
        μy = typeof(y) <: Real ? zero(T) : mean(y)
        μz = typeof(z) <: Real ? zero(T) : mean(z)
        origin = CartesianPoint{T}(μx, μy, μz)
        hX = typeof(x) <: Real ? x : (x[2] - x[1])/2
        hY = typeof(y) <: Real ? y : (y[2] - y[1])/2
        hZ = typeof(z) <: Real ? z : (z[2] - z[1])/2
        hX, hY, hZ
    end
    box = Box{T, ClosedPrimitive}(
        hX = hX, 
        hY = hY, 
        hZ = hZ, 
        origin = origin,
        rotation = rotation
    )
    transform(box, transformations)
end

function Dictionary(b::Box{T})::OrderedDict{String, Any} where {T}
    dict = OrderedDict{String,Any}()
    dict["widths"] = [2*b.hX, 2*b.hY, 2*b.hZ]
    if b.origin != zero(CartesianVector{T}) dict["origin"] = b.origin end
    if b.rotation != one(SMatrix{3,3,T,9}) dict["rotation"] = Dictionary(b.rotation) end
    OrderedDict{String,Any}("box" => dict)
end

function vertices(b::Box{T}) where {T}
    return (
        b.rotation * SVector{3,T}(-b.hX, -b.hY, -b.hZ) .+ b.origin,
        b.rotation * SVector{3,T}(+b.hX, -b.hY, -b.hZ) .+ b.origin,
        b.rotation * SVector{3,T}(+b.hX, +b.hY, -b.hZ) .+ b.origin,
        b.rotation * SVector{3,T}(-b.hX, +b.hY, -b.hZ) .+ b.origin,
        b.rotation * SVector{3,T}(-b.hX, -b.hY, +b.hZ) .+ b.origin,
        b.rotation * SVector{3,T}(+b.hX, -b.hY, +b.hZ) .+ b.origin,
        b.rotation * SVector{3,T}(+b.hX, +b.hY, +b.hZ) .+ b.origin,
        b.rotation * SVector{3,T}(-b.hX, +b.hY, +b.hZ) .+ b.origin,
    )
end

function sample(b::Box{T})::Vector{CartesianPoint{T}} where {T} 
    [vertices(b)...]
end

function surfaces(b::Box{T}) where {T}
    vs = vertices(b)
    return (
        Quadrangle{T}((vs[1], vs[2], vs[3], vs[4])),
        Quadrangle{T}((vs[5], vs[6], vs[2], vs[1])),
        Quadrangle{T}((vs[8], vs[7], vs[6], vs[5])),
        Quadrangle{T}((vs[6], vs[7], vs[3], vs[2])),
        Quadrangle{T}((vs[7], vs[8], vs[4], vs[3])),
        Quadrangle{T}((vs[8], vs[5], vs[1], vs[4])),
    )
end

