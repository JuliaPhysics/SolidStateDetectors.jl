"""
    struct Box{T, CO} <: AbstractVolumePrimitive{T}

Volume primitive describing a three-dimensional [Box](@ref) with its surfaces
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

See also [Constructive Solid Geometry (CSG)](@ref).
"""
struct Box{T,CO} <: AbstractVolumePrimitive{T,CO}
    hX::T
    hY::T
    hZ::T
    
    origin::CartesianPoint{T}
    rotation::SMatrix{3,3,T,9}
end

#Type conversion happens here
function Box{T}(CO,hX, hY, hZ, origin, rotation) where {T}
    _hX = _csg_convert_args(T, hX)
    _hY = _csg_convert_args(T, hY)
    _hZ = _csg_convert_args(T, hZ)
    Box{T,CO}(_hX, _hY, _hZ, origin, rotation)
end

#Type promotion happens here
function Box(CO, hX::TX, hY::TY, hZ::TZ, origin::PT, rotation::ROT) where {TX, TY, TZ, PT, ROT}
    eltypes = _csg_get_promoted_eltype.((TX, TY, TZ, PT, ROT))
    T = float(promote_type(eltypes...))
    Box{T}(CO,T(hX), T(hY), T(hZ), origin, rotation)
end

function Box(::Type{CO}=ClosedPrimitive;
    hX = 1,
    hY = 1,
    hZ = 1,
    origin = zero(CartesianPoint{Int}), 
    rotation = one(SMatrix{3, 3, Int, 9})
) where {CO}
    Box(CO, hX, hY, hZ, origin, rotation)
end

function Box{T}(::Type{CO}=ClosedPrimitive;
    hX = 1.0,
    hY = 1.0,
    hZ = 1.0,
    origin = zero(CartesianPoint{Float64}), 
    rotation = one(SMatrix{3, 3, Float64, 9})
) where {T, CO}
    Box{T}(CO, hX, hY, hZ, origin, rotation)
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
    end
    box = Box{T}(ClosedPrimitive,
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
    # dict["widths"] = [2*b.hX, 2*b.hY, 2*b.hZ]
    dict["hX"] = b.hX
    dict["hY"] = b.hY
    dict["hZ"] = b.hZ
    if !iszero(b.origin) dict["origin"] = Dictionary(b.origin) end
    if !isone(b.rotation) dict["rotation"] = Dictionary(b.rotation) end
    OrderedDict{String,Any}("box" => dict)
end

function vertices(b::Box{T}) where {T}
    return _transform_into_global_coordinate_system.((
        CartesianPoint{T}(-b.hX, -b.hY, -b.hZ),
        CartesianPoint{T}(+b.hX, -b.hY, -b.hZ),
        CartesianPoint{T}(+b.hX, +b.hY, -b.hZ),
        CartesianPoint{T}(-b.hX, +b.hY, -b.hZ),
        CartesianPoint{T}(-b.hX, -b.hY, +b.hZ),
        CartesianPoint{T}(+b.hX, -b.hY, +b.hZ),
        CartesianPoint{T}(+b.hX, +b.hY, +b.hZ),
        CartesianPoint{T}(-b.hX, +b.hY, +b.hZ)
    ), Ref(b))
end

function sample(b::Box{T})::Vector{CartesianPoint{T}} where {T} 
    [vertices(b)...]
end

function surfaces(b::Box{T,ClosedPrimitive}) where {T}
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

extremum(b::Box{T}) where {T} = hypot(b.hX, b.hY, b.hZ)