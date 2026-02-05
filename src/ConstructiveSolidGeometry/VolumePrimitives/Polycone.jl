"""
    struct Polyone{T,CO,N,TP} <: AbstractVolumePrimitive{T, CO}

Volume primitive describing a [Polycone](@ref), similar to the G4Polycone defined in Geant4.

## Parametric types
* `T`: Precision type.
* `CO`: Describes whether the surface belongs to the primitive. 
    It can be `ClosedPrimitive`, i.e. the surface points belong to the primitive,
    or `OpenPrimitive`, i.e. the surface points do not belong to the primitive.
* `N`: Integer describing the number of corners.
* `TP`: Type of the angular range `φ`.
    * `TP == Nothing`: Full 2π Cone.
    
## Fields
* `r::NTuple{N,T}`: `r`-coordinates of the corners of the polycone.
* `z::NTuple{N,T}`: `z`-coordinates of the corners of the polycone.
* `φ::TP`: Range in polar angle `φ` over which the `Cone` extends (in radians).
* `origin::CartesianPoint{T}`: The position of the center of the `Cone` at the middle height.
* `rotation::SMatrix{3,3,T,9}`: Matrix that describes a rotation of the `Cone` around its `origin`.

## Definition in Configuration File

A `Polycone` is defined in the configuration file as part of the `geometry` field 
of an object through the field `polycone`.

Example definitions of a polycone looks like this:
```yaml
polycone:
  r: [0, 35, 35, 24.42, 5, 5, 0, 0]
  z: [0, 0, 20, 80, 80, 25, 25, 0]
  origin:
    z: 1.0   # => origin = [0.0, 0.0, 1.0]
```
This is a polycone describing the semiconductor of the example inverted coaxial detector.

See also [Constructive Solid Geometry (CSG)](@ref).
"""
struct Polycone{T,CO,N,TP<:Nothing} <: AbstractVolumePrimitive{T,CO}
    r::NTuple{N,T}
    z::NTuple{N,T}
    φ::TP
    
    origin::CartesianPoint{T}
    rotation::SMatrix{3,3,T,9}
    
    function Polycone{T,CO}(r::Union{<:AbstractArray, <:Tuple}, z::Union{<:AbstractArray, <:Tuple}, φ, origin, rotation) where {T, CO}
        nr = length(r)
        nz = length(z)
        nr != nz && throw(ConfigFileError("In PolyCone: r and z must have the same length."))
        # convert to Tuple and determine the length N
        _r, _z = if first(r) != last(r) || first(z) != last(z)
            nr += 1
            nz += 1
            (r..., first(r)), (z..., first(z))
        else
            (r...,), (z...,)
        end
        # sort the points counter-clockwise in the r-z-plane
        if PolygonOps.area(tuple.(_r,_z)) < 0
            _r, _z = reverse(_r), reverse(_z)
        end
        new{T,CO,nr,typeof(φ)}(_r, _z, φ, origin, rotation)
    end
end

#Type conversion happens here
function Polycone{T}(CO, r, z, φ, origin, rotation) where {T}
    _r = _csg_convert_args.(T, r)
    _z = _csg_convert_args.(T, z)
    _φ = _csg_convert_args(T, φ)
    Polycone{T,CO}(_r, _z, _φ, origin, rotation)
end

#Type promotion happens here
function Polycone(CO, r::TR, z::TZ, φ::TP, origin::PT, rotation::ROT) where {TR, TZ, TP, PT, ROT}
    eltypes = _csg_get_promoted_eltype.((TR, TZ, TP, PT, ROT))
    T = float(promote_type(eltypes...))
    Polycone{T}(CO, r, z, φ, origin, rotation)
end

function Polycone(::Type{CO}=ClosedPrimitive;
    r = (0,1,1,0),
    z = (0,0,1,1),
    φ = nothing,
    origin = zero(CartesianPoint{Int}), 
    rotation = one(SMatrix{3, 3, Int, 9})
) where {CO}
    Polycone(CO, r, z, φ, origin, rotation)
end

function Polycone{T}(::Type{CO}=ClosedPrimitive;
    r = (0.,1.,1.,0.),
    z = (0.,0.,1.,1.),
    φ = nothing,
    origin = zero(CartesianPoint{Float64}), 
    rotation = one(SMatrix{3, 3, Float64, 9})
) where {T, CO}
    Polycone{T}(CO, r, z, φ, origin, rotation)
end

Polycone{T,CO,N,TP}( c::Polycone{T,CO,N,TP}; COT = CO,
            origin::CartesianPoint{T} = c.origin,
            rotation::SMatrix{3,3,T,9} = c.rotation) where {T, CO<:Union{ClosedPrimitive, OpenPrimitive}, N, TP} =
    Polycone{T,COT}(c.r, c.z, c.φ, origin, rotation)


####################################################################
####################################################################


function _inpolygon(pt::Tuple{T,T}, polygon::NTuple{N,Tuple{T,T}}; csgtol::T = csg_default_tol(T)) where {N,T <: Real}
    # use the normal PolygonOps.inpolygon method to check if the point is inside
    PolygonOps.inpolygon(pt, polygon, in = true, on = true, out = false) && return true
    # if it is flagged as outside, check if the distance to any of the edges is smaller than csgtol
    rp,zp = pt
    @inbounds for i in Base.OneTo(N-1)
        r1,z1 = polygon[i]
	r2,z2 = polygon[i+1]
        r1 == r2 && z1 == z2 && continue
        s = clamp(((zp - z1) * (z2 - z1) + (rp - r1) * (r2 - r1)) / ((z2 - z1)^2 + (r2 - r1)^2), zero(T), one(T))
        if hypot(rp - r1 - s * (r2 - r1), zp - z1 - s * (z2 - z1)) <= csgtol
            return true
        end
    end
    # if not, then it is really outside
    return false
end

function _in(pt::CartesianPoint{T}, c::Polycone{<:Any, ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T}
    return (isnothing(c.φ) || _in_angular_interval_closed(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)) && 
        _inpolygon((hypot(pt.x, pt.y), pt.z), tuple.(c.r, c.z); csgtol) 
end

function _in(pt::CartesianPoint{T}, c::Polycone{<:Any, OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T}
    return (isnothing(c.φ) || _in_angular_interval_open(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)) && 
        PolygonOps.inpolygon((hypot(pt.x, pt.y), pt.z), tuple.(c.r, c.z), in = true, on = iszero(pt.x) && iszero(pt.y), out = false)
end

function surfaces(c::Polycone{T,<:Any,N,Nothing}) where {T,N}
    s = []
    @inbounds for i in Base.OneTo(N-1)
        r1::T = c.r[i]
        r2::T = c.r[i+1]
        ## skip the surfaces on the z-axis
        iszero(r1) && iszero(r2) && continue
        z1::T = c.z[i]
        z2::T = c.z[i+1]
        origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), (z1+z2)/2), c) 
        vol = if z1 == z2
            EllipticalSurface{T}(r = r1 <= r2 ? (r1,r2) : (r2,r1), origin = origin, rotation = r1 <= r2 ? c.rotation : -c.rotation * RotZ{T}(π))
        else
            FullConeMantle{T, z1 < z2 ? (:inwards) : (:outwards)}(z1 < z2 ? (r1,r2) : (r2,r1), c.φ, abs(z1-z2)/2, origin, c.rotation)
        end
        push!(s, vol)
    end
    (s...,)
end


function Geometry(::Type{T}, ::Type{Polycone}, dict::AbstractDict, input_units::NamedTuple, transformations::Transformations{T}) where {T}
    length_unit = input_units.length
    origin = get_origin(T, dict, length_unit)
    rotation = get_rotation(T, dict, input_units.angle)
    haskey(dict, "r") || throw(ConfigFileError("Polycone needs entry `r`."))
    haskey(dict, "z") || throw(ConfigFileError("Polycone needs entry `z`."))
    r = _parse_value(T, dict["r"], length_unit)
    z = _parse_value(T, dict["z"], length_unit)
    polycone = Polycone{T, ClosedPrimitive}(r, z, nothing, origin, rotation)
    transform(polycone, transformations)
end


function Dictionary(c::Polycone{T})::OrderedDict{String, Any} where {T}
    dict = OrderedDict{String, Any}()
    @assert isnothing(c.φ) "Polycone needs `φ` field to be nothing."
    dict["r"] = c.r
    dict["z"] = c.z
    if !iszero(c.origin) dict["origin"] = Dictionary(c.origin) end
    if !isone(c.rotation) dict["rotation"] = Dictionary(c.rotation) end 
    OrderedDict{String, Any}("polycone" => dict)
end

extremum(c::Polycone{T}) where {T} = hypot(max(abs.(extrema(c.r))...), max(abs.(extrema(c.z))...))
