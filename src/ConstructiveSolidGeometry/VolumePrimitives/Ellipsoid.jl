"""
    struct Ellipsoid{T,CO,TR,TP,TT} <: AbstractVolumePrimitive{T, CO}

Volume primitive describing an [Ellipsoid](@ref).

## Parametric types
* `T`: Precision type.
* `CO`: Describes whether the surface belongs to the primitive. 
    It can be `ClosedPrimitive`, i.e. the surface points belong to the primitive,
    or `OpenPrimitive`, i.e. the surface points do not belong to the primitive.
* `TR`: Type of the radius `r`.
    * `TR == T`: Sphere (constant radius `r` along all axes).
* `TP`: Type of the azimuthial angle `φ`.
    * `TP == Nothing`: Full 2π in `φ`.
* `TT`: Type of the polar angle `θ`.
    * `TT == Nothing`: Full 2π in `θ`.

## Fields
* `r::TR`: Definition of the radius of the `Ellipsoid` (in m).
* `φ::TP`: Range in azimuthial angle `φ` of the `Ellipsoid`.
* `θ::TT`: Range in polar angle `θ` of the `Ellipsoid`.
* `origin::CartesianPoint{T}`: The position of the center of the `Ellipsoid`.
* `rotation::SMatrix{3,3,T,9}`: Matrix that describes a rotation of the `Ellipsoid` around its `origin`.

## Definition in Configuration File

So far, the only `Ellipsoid` implemented so far is a `FullSphere`.
A `FullSphere` is defined in the configuration file as part of the `geometry` field 
of an object through the field `sphere`.

Example definitions of a `FullSphere` looks like this:
```yaml
sphere:
  r: 2
```
This is a full sphere with radius 2.

To define a sphere with inner cut-out, use [`CSGDifference`](@ref):
```yaml
difference:
  - sphere:
      r: 2
  - sphere:
      r: 1
```
This is a sphere with inner radius 1 and outer radius 2.

See also [Constructive Solid Geometry (CSG)](@ref).
"""
struct Ellipsoid{T,CO,TR,TP,TT} <: AbstractVolumePrimitive{T, CO}
    r::TR 
    φ::TP 
    θ::TT 

    origin::CartesianPoint{T}
    rotation::SMatrix{3,3,T,9}    
end

#Type conversion happens here
function Ellipsoid{T,CO}(r, origin, rotation) where {T,CO}
    _r = _csg_convert_args(T, r)
    Ellipsoid{T,CO,typeof(_r),Nothing,Nothing}(_r, nothing, nothing, origin, rotation)
end

#Type promotion happens here
function Ellipsoid(CO, r::TR, origin::PT, rotation::ROT) where {TR, PT, ROT}
    eltypes = _csg_get_promoted_eltype.((TR, PT, ROT))
    T = float(promote_type(eltypes...))
    Ellipsoid{T,CO}(r, origin, rotation)
end

function Ellipsoid(::Type{CO}=ClosedPrimitive;
    r = 1, 
    origin = zero(CartesianPoint{Int}), 
    rotation = one(SMatrix{3, 3, Int, 9})
) where {CO}
    Ellipsoid(CO, r, origin, rotation)
end

function Ellipsoid{T}(::Type{CO}=ClosedPrimitive;
    r = 1.0, 
    origin = zero(CartesianPoint{Float64}), 
    rotation = one(SMatrix{3, 3, Float64, 9})
) where {T, CO}
    Ellipsoid{T,CO}( r, origin, rotation)
end

Ellipsoid{T,CO,TR,TP,TT}( e::Ellipsoid{T,CO,TR,TP,TT}; COT = CO,
            origin::CartesianPoint{T} = e.origin,
            rotation::SMatrix{3,3,T,9} = e.rotation) where {T,CO<:Union{ClosedPrimitive, OpenPrimitive},TR,TP,TT} =
    Ellipsoid{T,COT,TR,TP,TT}(e.r,e.φ,e.θ,origin,rotation)

const Sphere{T,CO,TP,TT} = Ellipsoid{T,CO,T,TP,TT}
const FullSphere{T,CO} = Ellipsoid{T,CO,T,Nothing,Nothing}
const FullEllipsoid{T,CO} = Ellipsoid{T,CO,NTuple{3,T},Nothing,Nothing}

function Geometry(::Type{T}, ::Type{Ellipsoid}, dict::AbstractDict, input_units::NamedTuple, transformations::Transformations{T}) where {T}
    length_unit = input_units.length
    angle_unit = input_units.angle
    origin = get_origin(T, dict, length_unit)
    rotation = get_rotation(T, dict, angle_unit)
    
    r = parse_r_of_primitive(T, dict, input_units.length)
    φ = parse_φ_of_primitive(T, dict, angle_unit)
    if !(φ === nothing)
        error("Partial Ellipsoid (`φ = φ`) is not yet supported.")
    end
    θ = parse_θ_of_primitive(T, dict, angle_unit)
    if !(θ === nothing)
        error("Partial Ellipsoid (`θ = θ`) is not yet supported.")
    end
    e = if r isa Tuple{T,T}
        Ellipsoid{T}(ClosedPrimitive,
            r = r[2], 
            origin = origin,
            rotation = rotation,
        ) - Ellipsoid{T}(OpenPrimitive,
            r = r[1], 
            origin = origin,
            rotation = rotation,
        )
    else
        Ellipsoid{T}(ClosedPrimitive,
            r = r, 
            origin = origin,
            rotation = rotation,
        )
    end
    transform(e, transformations)
end

function Dictionary(e::Ellipsoid{T})::OrderedDict{String, Any} where {T}
    dict = OrderedDict{String, Any}()
    dict["r"] = e.r # always a Real 
    if !isnothing(e.φ) error("Partial Ellipsoid (`φ = φ`) is not yet supported.") end
    if !isnothing(e.θ) error("Partial Ellipsoid (`θ = θ`) is not yet supported.") end
    if e.origin != zero(CartesianVector{T}) dict["origin"] = e.origin end
    if e.rotation != one(SMatrix{3,3,T,9}) dict["rotation"] = Dictionary(e.rotation) end
    OrderedDict{String, Any}("sphere" => dict)
end

_in(pt::CartesianPoint{T}, s::FullSphere{<:Any, ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} =
    hypot(pt.x, pt.y, pt.z) <= s.r + csgtol

_in(pt::CartesianPoint{T}, s::FullSphere{<:Any, OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} =
    hypot(pt.x, pt.y, pt.z) < s.r - csgtol

function surfaces(e::Ellipsoid{T,ClosedPrimitive}) where {T}
    em = EllipsoidMantle{T,typeof(e.r),typeof(e.φ),typeof(e.θ),:inwards}(e.r, e.φ, e.θ, e.origin, e.rotation)
    (em,)
end

extremum(e::Ellipsoid{T}) where {T} = max(e.r...)