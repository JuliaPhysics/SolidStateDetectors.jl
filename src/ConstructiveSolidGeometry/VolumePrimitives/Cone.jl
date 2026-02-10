"""
    struct Cone{T,CO,TR,TP} <: AbstractVolumePrimitive{T, CO}

Volume primitive describing a [Cone](@ref) with its top and bottom circular base
being aligned with the `xy` plane (before possible rotations).

## Parametric types
* `T`: Precision type.
* `CO`: Describes whether the surface belongs to the primitive. 
    It can be `ClosedPrimitive`, i.e. the surface points belong to the primitive,
    or `OpenPrimitive`, i.e. the surface points do not belong to the primitive.
* `TR`: Type of the radius `r`.
    * `TR == T`: Cylinder (constant radius `r` at all `z`).
    * `TR == Tuple{T, T}`: Tube (inner radius at `r[1]`, outer radius at `r[2]`).
    * `TR == Tuple{Tuple{T}, Tuple{T}}`: Varying Cylinder (full cylinder with radius changing linearly in `z` from `r[1][1]` at the bottom to `r[2][1]` at the top).
    * `TR == Tuple{Tuple{T, T}, Tuple{T, T}}`: Varying Tube (inner radius changes linearly in `z` from `r[1][1]` at the bottom to `r[2][1]` at the top, outer radius changes linearly in `z` from `r[1][2]` at the bottom to `r[2][2]` at the top).
    * `TR == Tuple{Nothing, Tuple{T, T}}`: Cone (Tip at the bottom, top is a circular base with inner radius `r[2][1]` and outer radius `r[2][2]`).
    * `TR == Tuple{Tuple{T, T}, Nothing}`: Cone (Tip at the top, bottom is a circular base with inner radius `r[1][1]` and outer radius `r[1][2]`).
* `TP`: Type of the angular range `φ`.
    * `TP == Nothing`: Full 2π Cone.
    * `TP == T`: Partial Cone ranging from `0` to `φ`.
    
## Fields
* `r::TR`: Definition of the radius of the `Cone` (in m).
* `φ::TP`: Range in polar angle `φ` over which the `Cone` extends (in radians).
* `hZ::T`: Half of the height of the `Cone` (in m).
* `origin::CartesianPoint{T}`: The position of the center of the `Cone` at the middle height.
* `rotation::SMatrix{3,3,T,9}`: Matrix that describes a rotation of the `Cone` around its `origin`.

## Definition in Configuration File

A `Cone` is defined in the configuration file as part of the `geometry` field 
of an object through the field `cone` (or `tube`).

Example definitions of a cylinder looks like this:
```yaml
tube:
  r:
    from: 1.0
    to: 2.0  # => r = (1.0, 2.0)
  h: 2.0     # => hZ = 1.0
  origin:
    z: 1.0   # => origin = [0.0, 0.0, 1.0]
```
This is a hollow cylinder with inner radius 1 at outer radius 2
with a height of 2 and extending over full 2π (no `phi` given),
translated 1 along the `z` axis.

If the radius is not constant over `z`, the `r` entries are divided
into `bottom` and `top`, where `bottom` describes the inner and outer 
radius at the bottom circular base and `top` describes the inner and
outer radius at the top circular base (before rotations)
```yaml
cone:
  r:
    bottom:
      from: 1.0
      to: 2.0
    top:
      from: 1.0
      to: 4.0     # => r = ((1.0, 2.0), (1.0, 4.0))
  phi:
    from: 0.0°
    to: 180.0°    # => φ = π
  h: 2.0          # => hZ = 1.0
```
This is half a Cone (`φ` goes from 0 to 180°, i.e. only positive `y` are allowed)
with a height of 2, constant inner radius of 1 and an outer radius which increases 
from 2 at the bottom to 4 at the top circular base.

!!! note
    The names `tube` and `cone` in the configuration files are interchangeable.
    
See also [Constructive Solid Geometry (CSG)](@ref).
"""
struct Cone{T,CO,TR,TP<:Union{Nothing,T}} <: AbstractVolumePrimitive{T, CO}
    r::TR 
    φ::TP
    hZ::T 

    origin::CartesianPoint{T} 
    rotation::SMatrix{3,3,T,9}
end

#Type conversion happens here
function Cone{T,CO}(r, φ, hZ, origin, rotation) where {T,CO}
    _r = _csg_convert_args(T, r)
    (_φ, _rotation) = _handle_phi(_csg_convert_args(T, φ), rotation)
    _hZ = _csg_convert_args(T, hZ)
    Cone{T,CO,typeof(_r),typeof(_φ)}(_r, _φ, _hZ, origin, _rotation)
end

#Type promotion happens here
function Cone(CO, r::TR, φ::TP, hZ::TZ, origin::PT, rotation::ROT) where {TR, TP, TZ, PT, ROT}
    eltypes = _csg_get_promoted_eltype.((TR, TZ, TP, PT, ROT))
    T = float(promote_type(eltypes...))
    Cone{T,CO}(r, φ, hZ, origin, rotation)
end

function Cone(::Type{CO}=ClosedPrimitive;
    # define default parameters as Int to not influence type promotion
    r = 1, 
    φ = nothing,
    hZ = 1,
    origin = zero(CartesianPoint{Int}), 
    rotation = one(SMatrix{3, 3, Int, 9})
) where {CO}
    Cone(CO, r, φ, hZ, origin, rotation)
end

function Cone{T}(::Type{CO}=ClosedPrimitive;
    r = 1.0, 
    φ = nothing,
    hZ = 1.0,
    origin = zero(CartesianPoint{Float64}), 
    rotation = one(SMatrix{3, 3, Float64, 9})
) where {T, CO}
    Cone{T,CO}(r, φ, hZ, origin, rotation)
end

Cone{T,CO,TR,TP}( c::Cone{T,CO,TR,TP}; COT = CO,
            origin::CartesianPoint{T} = c.origin,
            rotation::SMatrix{3,3,T,9} = c.rotation) where {T,CO<:Union{ClosedPrimitive, OpenPrimitive},TR,TP<:Union{Nothing,T}} =
    Cone{T,COT,TR,TP}(c.r, c.φ, c.hZ, origin, rotation)

####################################################################
####################################################################


#=
___
   |
   |
   |
___|

=#
const Cylinder{T,CO} = Cone{T,CO,T,Nothing} # Full in φ
const PartialCylinder{T,CO} = Cone{T,CO,T,T} 

### Cylinder
function _in(pt::CartesianPoint, c::Cylinder{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    hypot(max(zero(T), hypot(pt.x, pt.y) - c.r), max(zero(T), pt.z - c.hZ, - pt.z - c.hZ)) <= csgtol
end

function _in(pt::CartesianPoint, c::Cylinder{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) < c.hZ - csgtol && hypot(pt.x, pt.y) < c.r - csgtol
end

function surfaces(t::Cylinder{T,ClosedPrimitive}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    mantle = FullConeMantle{T,:inwards}((t.r, t.r), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = EllipticalSurface{T}(r = t.r, φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = EllipticalSurface{T}(r = t.r, φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    # normals of the surfaces show inside the volume primitives. 
    e_top, e_bot, mantle
end

### PartialCylinder
function _in(pt::CartesianPoint, c::PartialCylinder{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T}
    r::T = hypot(pt.x, pt.y)
    φ::T = mod(atan(pt.y, pt.x), T(2π))
    z::T = pt.z
    abs(z) <= c.hZ + csgtol && begin
        det::T = csgtol^2 - max(zero(T), z - c.hZ, - z - c.hZ)^2
        det >= 0 && (_in_angular_interval_closed(φ, c.φ, csgtol = zero(T)) && r <= c.r + sqrt(det) || let Δφ = min(T(2π)-φ, φ-c.φ)
            if Δφ <= atan(sqrt(det), c.r)
                det >= c.r^2 * sin(Δφ)^2 && r <= c.r * cos(Δφ) + sqrt(det - c.r^2 * sin(Δφ)^2)
            else
                r <= sqrt(det)/sin(min(Δφ, T(π/2)))
            end
        end)
    end
end

function _in(pt::CartesianPoint, c::PartialCylinder{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T}
    r::T = hypot(pt.x, pt.y)
    φ::T = mod(atan(pt.y, pt.x), T(2π))
    abs(pt.z) + csgtol < c.hZ && r < c.r - csgtol && _in_angular_interval_closed(φ, c.φ, csgtol = zero(T)) && r * sin(min(φ, c.φ - φ, T(π/2))) > csgtol
end

function surfaces(t::PartialCylinder{T,ClosedPrimitive}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    mantle = PartialConeMantle{T,:inwards}((t.r, t.r), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = EllipticalSurface{T}(r = t.r, φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = EllipticalSurface{T}(r = t.r, φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    poly_l = _transform_into_global_coordinate_system(Quadrangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T), zero(T), -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T), zero(T), +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(    t.r, zero(T), +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(    t.r, zero(T), -t.hZ)) )), t)
    poly_r = _transform_into_global_coordinate_system(Quadrangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T), t.φ, -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(    t.r, t.φ, -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(    t.r, t.φ, +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T), t.φ, +t.hZ)) )), t)
    # normals of the surfaces show inside the volume primitives. 
    (e_top, e_bot, mantle, poly_l, poly_r)
end


####################################################################
####################################################################

#=
___
   \
    \
_____\

=#
const VaryingCylinder{T,CO} = Cone{T,CO,Tuple{Tuple{T},Tuple{T}},Nothing} # Full in φ
const PartialVaryingCylinder{T,CO} = Cone{T,CO,Tuple{Tuple{T},Tuple{T}},T}


### VaryingCylinder
function _in(pt::CartesianPoint, c::VaryingCylinder{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    r::T = hypot(pt.x, pt.y)
    z::T = pt.z 
    rz::T = radius_at_z(c.hZ, c.r[1][1], c.r[2][1], z)
    Δr::T = c.r[2][1] - c.r[1][1]
    Δz::T = iszero(c.hZ) ? zero(T) : -csgtol * Δr / hypot(2*c.hZ, Δr)
    return (abs(z - Δz) < c.hZ && r <= rz + csgtol * hypot(2*c.hZ, Δr) / (2*c.hZ)) ||
        ((z + c.hZ)^2 <= csgtol^2 && r <= c.r[1][1] + sqrt(csgtol^2 - (z + c.hZ)^2)) ||
        ((z - c.hZ)^2 <= csgtol^2 && r <= c.r[2][1] + sqrt(csgtol^2 - (z - c.hZ)^2))
end

function _in(pt::CartesianPoint, c::VaryingCylinder{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) < c.hZ - csgtol &&
    hypot(pt.x, pt.y) < radius_at_z(c.hZ, c.r[1][1], c.r[2][1], pt.z) - csgtol * hypot(2*c.hZ, c.r[2][1] - c.r[1][1]) / (2*c.hZ)
end

function surfaces(t::VaryingCylinder{T,ClosedPrimitive}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    mantle = FullConeMantle{T,:inwards}((t.r[1][1], t.r[2][1]), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = EllipticalSurface{T}(r = t.r[1][1], φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = EllipticalSurface{T}(r = t.r[2][1], φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    # normals of the surfaces show inside the volume primitives. 
    e_top, e_bot, mantle
end

### PartialVaryingCylinder
function _in(pt::CartesianPoint, c::PartialVaryingCylinder{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    r::T = hypot(pt.x, pt.y)
    φ::T = mod(atan(pt.y, pt.x), T(2π))
    z::T = pt.z
    rz::T = radius_at_z(c.hZ, c.r[1][1], c.r[2][1], z)
    Δr::T = c.r[2][1] - c.r[1][1]
    abs(z) <= c.hZ + csgtol && if _in_angular_interval_closed(φ, c.φ, csgtol = zero(T)) 
        _in(pt, Cone{T,ClosedPrimitive}(c.r, nothing, c.hZ, c.origin, c.rotation); csgtol)
    else
        t::T = ((r - c.r[1][1]) * Δr + (z + c.hZ) * 2c.hZ) / (Δr^2 + 4c.hZ^2)
        s::T = ((r - c.r[1][1]) * 2c.hZ - (z + c.hZ) * Δr) / (Δr^2 + 4c.hZ^2)
        d::T = if r <= rz && abs(z) <= c.hZ
            zero(T)
        elseif 0 <= t <= 1 && s >= 0
            s * hypot(Δr, 2c.hZ)
        elseif t > 1 || (t >= 0 && s <= 0 && Δr > 0)
            hypot(abs(z - c.hZ), max(zero(T), r - c.r[2][1]))
        elseif t < 0 || (t <= 1 && s <= 0 && Δr < 0)
            hypot(abs(-z - c.hZ), max(zero(T), r - c.r[1][1]))
        else 
            hypot(r, max(zero(T), pt.z - c.hZ, - pt.z - c.hZ))
        end
        (r * sin(min(T(2π)-φ, φ-c.φ, T(π/2))))^2 + d^2 <= csgtol^2
    end
end

function _in(pt::CartesianPoint, c::PartialVaryingCylinder{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    r::T = hypot(pt.x, pt.y) 
    φ::T = mod(atan(pt.y, pt.x), T(2π))
    Δr::T = c.r[2][1] - c.r[1][1]
    abs(pt.z) + csgtol < c.hZ && 
        r < radius_at_z(c.hZ, c.r[1][1], c.r[2][1], pt.z) - csgtol * hypot(2c.hZ, Δr) / 2c.hZ &&
        _in_angular_interval_open(atan(pt.y, pt.x), c.φ, csgtol = zero(T)) && r * sin(min(φ, c.φ - φ, T(π/2))) > csgtol
end

function surfaces(t::PartialVaryingCylinder{T,ClosedPrimitive}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    mantle = PartialConeMantle{T,:inwards}((t.r[1][1], t.r[2][1]), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = EllipticalSurface{T}(r = t.r[1][1], φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = EllipticalSurface{T}(r = t.r[2][1], φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    poly_l = _transform_into_global_coordinate_system(Quadrangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T),    zero(T), -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T),    zero(T), +t.hZ)),
        CartesianPoint(CylindricalPoint{T}( t.r[2][1], zero(T), +t.hZ)),
        CartesianPoint(CylindricalPoint{T}( t.r[1][1], zero(T), -t.hZ)) )), t)
    poly_r = _transform_into_global_coordinate_system(Quadrangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T),    t.φ, -t.hZ)),
        CartesianPoint(CylindricalPoint{T}( t.r[1][1], t.φ, -t.hZ)),
        CartesianPoint(CylindricalPoint{T}( t.r[2][1], t.φ, +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T),    t.φ, +t.hZ)) )), t)
    # normals of the surfaces show inside the volume primitives. 
    e_top, e_bot, mantle, poly_l, poly_r
end



####################################################################
####################################################################

#=
  ___
 /   \
/_____\

=#
const VaryingTube{T,CO} = Cone{T,CO,Tuple{Tuple{T,T},Tuple{T,T}},Nothing} # Full in φ
const PartialVaryingTube{T,CO} = Cone{T,CO,Tuple{Tuple{T,T},Tuple{T,T}},T}      

#(r_bot_in = r[1][1], r_bot_out = r[1][2], r_top_in = r[2][1], r_top_out = r[2][2])

### VaryingTube
function _in(pt::CartesianPoint, c::VaryingTube{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T}
    r::T = hypot(pt.x, pt.y)
    z::T = pt.z 
    r_in::T  = radius_at_z(c.hZ, c.r[1][1], c.r[2][1], z)
    r_out::T = radius_at_z(c.hZ, c.r[1][2], c.r[2][2], z)
    Δr_in::T  = c.r[2][1] - c.r[1][1]
    Δr_out::T = c.r[2][2] - c.r[1][2]
    Δz_in::T  = iszero(c.hZ) ? zero(T) : csgtol * Δr_in / hypot(2*c.hZ, Δr_in)
    Δz_out::T = iszero(c.hZ) ? zero(T) : csgtol * Δr_out / hypot(2*c.hZ, Δr_out)
    return (abs(z - Δz_in) < c.hZ && r >= r_in - csgtol * hypot(2*c.hZ, Δr_in) / (2*c.hZ) ||
        ((z + c.hZ)^2 <= csgtol^2 && r >= c.r[1][1] - sqrt(csgtol^2 - (z + c.hZ)^2)) ||
        ((z - c.hZ)^2 <= csgtol^2 && r >= c.r[2][1] - sqrt(csgtol^2 - (z - c.hZ)^2))) &&
        (abs(z + Δz_out) < c.hZ && r <= r_out + csgtol * hypot(2*c.hZ, Δr_out) / (2*c.hZ) ||
        ((z + c.hZ)^2 <= csgtol^2 && r <= c.r[1][2] + sqrt(csgtol^2 - (z + c.hZ)^2)) ||
        ((z - c.hZ)^2 <= csgtol^2 && r <= c.r[2][2] + sqrt(csgtol^2 - (z - c.hZ)^2)))
end

function _in(pt::CartesianPoint, c::VaryingTube{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T}
    r::T = hypot(pt.x, pt.y)
    return abs(pt.z) + csgtol < c.hZ &&
        r > radius_at_z(c.hZ, c.r[1][1], c.r[2][1], pt.z) + csgtol * hypot(2*c.hZ, c.r[2][1] - c.r[1][1]) / (2*c.hZ) &&
        r < radius_at_z(c.hZ, c.r[1][2], c.r[2][2], pt.z) - csgtol * hypot(2*c.hZ, c.r[2][2] - c.r[1][2]) / (2*c.hZ)
end

function surfaces(t::VaryingTube{T,ClosedPrimitive}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    inner_mantle = FullConeMantle{T,:outwards}((t.r[1][1], t.r[2][1]), t.φ, t.hZ, t.origin, t.rotation)
    outer_mantle = FullConeMantle{T,:inwards}( (t.r[1][2], t.r[2][2]), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = EllipticalSurface{T}(r = (t.r[1][1], t.r[1][2]), φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = EllipticalSurface{T}(r = (t.r[2][1], t.r[2][2]), φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    # normals of the surfaces show inside the volume primitives. 
    e_top, e_bot, inner_mantle, outer_mantle
end

### PartialVaryingTube
function _in(pt::CartesianPoint, c::PartialVaryingTube{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    r::T = hypot(pt.x, pt.y)
    φ::T = mod(atan(pt.y, pt.x), T(2π))
    z::T = pt.z
    abs(z) <= c.hZ + csgtol && if _in_angular_interval_closed(atan(pt.y, pt.x), c.φ, csgtol = zero(T))
        _in(pt, Cone{T,ClosedPrimitive}(c.r, nothing, c.hZ, c.origin, c.rotation); csgtol)
    else
        r_in  = clamp(radius_at_z(c.hZ, c.r[1][1], c.r[2][1], z), extrema((c.r[1][1], c.r[2][1]))...)
        r_out = clamp(radius_at_z(c.hZ, c.r[1][2], c.r[2][2], z), extrema((c.r[1][2], c.r[2][2]))...)
        Δr_in  = c.r[2][1] - c.r[1][1]
        Δr_out = c.r[2][2] - c.r[1][2]
        t_in =  ((r - c.r[1][1]) * Δr_in + (z + c.hZ) * 2c.hZ) / (Δr_in^2 + 4c.hZ^2)
        s_in = -((r - c.r[1][1]) * 2c.hZ - (z + c.hZ) * Δr_in) / (Δr_in^2 + 4c.hZ^2)
        t_out = ((r - c.r[1][2]) * Δr_out + (z + c.hZ) * 2c.hZ) / (Δr_out^2 + 4c.hZ^2)
        s_out = ((r - c.r[1][2]) * 2c.hZ - (z + c.hZ) * Δr_out) / (Δr_out^2 + 4c.hZ^2)
        d::T = if r_in <= r <= r_out && abs(z) <= c.hZ
            zero(T)
        elseif 0 <= t_in <= 1 && s_in >= 0
            s_in * hypot(Δr_in, 2c.hZ)
        elseif 0 <= t_out <= 1 && s_out >= 0
            s_out * hypot(Δr_out, 2c.hZ)
        elseif t_out > 1 || t_in > 1
            hypot(abs(z - c.hZ), max(0, r - c.r[2][2], c.r[2][1] - r))
        else # t_out < 0 || t_in < 0
            hypot(abs(-z - c.hZ), max(0, r - c.r[1][2], c.r[1][1] - r))
        end
        (r * sin(min(T(2π)-φ, φ-c.φ, T(π/2))))^2 + d^2 <= csgtol^2
    end 
end

function _in(pt::CartesianPoint, c::PartialVaryingTube{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    r::T = hypot(pt.x, pt.y) 
    φ::T = mod(atan(pt.y, pt.x), T(2π))
    abs(pt.z) + csgtol < c.hZ && 
        _in_angular_interval_open(atan(pt.y, pt.x), c.φ, csgtol = zero(T)) && r * sin(min(φ, c.φ - φ, T(π/2))) > csgtol &&
        r > radius_at_z(c.hZ, c.r[1][1], c.r[2][1], pt.z) + csgtol * hypot(2c.hZ, c.r[2][1] - c.r[1][1]) / 2c.hZ &&
        r < radius_at_z(c.hZ, c.r[1][2], c.r[2][2], pt.z) - csgtol * hypot(2c.hZ, c.r[2][2] - c.r[1][2]) / 2c.hZ
end

function surfaces(t::PartialVaryingTube{T,ClosedPrimitive}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    inner_mantle = PartialConeMantle{T,:outwards}((t.r[1][1], t.r[2][1]), t.φ, t.hZ, t.origin, t.rotation)
    outer_mantle = PartialConeMantle{T,:inwards}( (t.r[1][2], t.r[2][2]), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = EllipticalSurface{T}(r = (t.r[1][1], t.r[1][2]), φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = EllipticalSurface{T}(r = (t.r[2][1], t.r[2][2]), φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    poly_l = _transform_into_global_coordinate_system(Quadrangle{T}((
        CartesianPoint(CylindricalPoint{T}(t.r[1][1], zero(T), -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][1], zero(T), +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][2], zero(T), +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[1][2], zero(T), -t.hZ)) )), t)
    poly_r = _transform_into_global_coordinate_system(Quadrangle{T}((
        CartesianPoint(CylindricalPoint{T}(t.r[1][1], t.φ, -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[1][2], t.φ, -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][2], t.φ, +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][1], t.φ, +t.hZ)) )), t)
    # normals of the surfaces show inside the volume primitives. 
    e_top, e_bot, inner_mantle, outer_mantle, poly_l, poly_r
end

####################################################################
####################################################################

#=

\
 \
__\

=#
const UpwardCone{T,CO} = Cone{T,CO,Tuple{Tuple{Nothing,T},Nothing},Nothing} # Full in φ
const PartialUpwardCone{T,CO} = Cone{T,CO,Tuple{Tuple{Nothing,T},Nothing},T}

function _in(pt::CartesianPoint, c::Cone{T,CO,Tuple{Tuple{Nothing,T},Nothing}}; csgtol::T = csg_default_tol(T)) where {T,CO} 
    _in(pt, Cone{T,CO}(((c.r[1][2],), (zero(T),)), c.φ, c.hZ, c.origin, c.rotation); csgtol)
end

function surfaces(t::UpwardCone{T,ClosedPrimitive}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    outer_mantle = FullConeMantle{T,:inwards}( (t.r[1][2], zero(T)), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = EllipticalSurface{T}(r = t.r[1][2], φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    # normals of the surfaces show inside the volume primitives. 
    e_bot, outer_mantle
end

function surfaces(t::PartialUpwardCone{T,ClosedPrimitive}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    outer_mantle = PartialConeMantle{T,:inwards}( (t.r[1][2], zero(T)), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = EllipticalSurface{T}(r = t.r[1][2], φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    poly_l = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T),   zero(T), -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T),   zero(T), +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[1][2], zero(T), -t.hZ)) )), t)
    poly_r = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T),   t.φ, -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[1][2], t.φ, -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T),   t.φ, +t.hZ)) )), t)
    # normals of the surfaces show inside the volume primitives. 
    e_bot, outer_mantle, poly_l, poly_r
end


####################################################################
####################################################################

#=
___
  /
 /
/

=#
const DownwardCone{T,CO} = Cone{T,CO,Tuple{Nothing,Tuple{Nothing,T}},Nothing} # Full in φ
const PartialDownwardCone{T,CO} = Cone{T,CO,Tuple{Nothing,Tuple{Nothing,T}},T} 

function _in(pt::CartesianPoint, c::Cone{T,CO,Tuple{Nothing,Tuple{Nothing,T}}}; csgtol::T = csg_default_tol(T)) where {T,CO} 
    _in(pt, Cone{T,CO}(((zero(T),), (c.r[2][2],)), c.φ, c.hZ, c.origin, c.rotation); csgtol)
end

function surfaces(t::DownwardCone{T,ClosedPrimitive}) where {T}
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    outer_mantle = FullConeMantle{T,:inwards}( (zero(T), t.r[2][2]), t.φ, t.hZ, t.origin, t.rotation)
    e_top = EllipticalSurface{T}(r = t.r[2][2], φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    # normals of the surfaces show inside the volume primitives. 
    e_top, outer_mantle
end

function surfaces(t::PartialDownwardCone{T,ClosedPrimitive}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    outer_mantle = PartialConeMantle{T,:inwards}( (t.r[2][2], zero(T)), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = EllipticalSurface{T}(r = t.r[2][2], φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    poly_l = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T),   zero(T), -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T),   zero(T), +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][2], zero(T), -t.hZ)) )), t)
    poly_r = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T),   t.φ, -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][2], t.φ, -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T),   t.φ, +t.hZ)) )), t)
    # normals of the surfaces show inside the volume primitives. 
    e_bot, outer_mantle, poly_l, poly_r
end



####################################################################
####################################################################

#=

  /\
 /  \
/____\

=#
const TopClosedTube{T,CO} = Cone{T,CO,Tuple{Tuple{T,T},T},Nothing} # Full in φ
const PartialTopClosedTube{T,CO} = Cone{T,CO,Tuple{Tuple{T,T},T},T} 

function _in(pt::CartesianPoint, c::Cone{T,CO,Tuple{Tuple{T,T},T}}; csgtol::T = csg_default_tol(T)) where {T,CO} 
    _in(pt, Cone{T,CO}((c.r[1], (c.r[2], c.r[2])), c.φ, c.hZ, c.origin, c.rotation); csgtol)
end

function surfaces(t::TopClosedTube{T,ClosedPrimitive}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    # top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    inner_mantle = FullConeMantle{T,:outwards}((t.r[1][1], t.r[2]), t.φ, t.hZ, t.origin, t.rotation)
    outer_mantle = FullConeMantle{T,:inwards}( (t.r[1][2], t.r[2]), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = EllipticalSurface{T}(r = (t.r[1][1], t.r[1][2]), φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    # e_top = Annulus{T}(r = (t.r[2][1], t.r[2][2]), φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    # normals of the surfaces show inside the volume primitives. 
    e_bot, inner_mantle, outer_mantle
end

function surfaces(t::PartialTopClosedTube{T,ClosedPrimitive}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    # top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    inner_mantle = PartialConeMantle{T,:outwards}((t.r[1][1], t.r[2]), t.φ, t.hZ, t.origin, t.rotation)
    outer_mantle = PartialConeMantle{T,:inwards}( (t.r[1][2], t.r[2]), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = EllipticalSurface{T}(r = (t.r[1][1], t.r[1][2]), φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    # e_top = Annulus{T}(r = (t.r[2][1], t.r[2][2]), φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    poly_l = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(t.r[1][1], zero(T), -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(   t.r[2], zero(T), +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[1][2], zero(T), -t.hZ)) )), t)
    poly_r = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(t.r[1][1], t.φ, -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[1][2], t.φ, -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(   t.r[2], t.φ, +t.hZ)) )), t)
    # normals of the surfaces show inside the volume primitives. 
    e_bot, inner_mantle, outer_mantle, poly_l, poly_r
end



####################################################################
####################################################################

#=

______
\    /
 \  /
  \/

=#
const BottomClosedTube{T,CO} = Cone{T,CO,Tuple{T,Tuple{T,T}},Nothing} # Full in φ
const PartialBottomClosedTube{T,CO} = Cone{T,CO,Tuple{T,Tuple{T,T}},T} 

function _in(pt::CartesianPoint, c::Cone{T,CO,Tuple{T,Tuple{T,T}}}; csgtol::T = csg_default_tol(T)) where {T,CO} 
    _in(pt, Cone{T,CO}(((c.r[1], c.r[1]), c.r[2]), c.φ, c.hZ, c.origin, c.rotation); csgtol)
end

function surfaces(t::BottomClosedTube{T,ClosedPrimitive}) where {T}
    # bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    inner_mantle = FullConeMantle{T,:outwards}((t.r[1], t.r[2][1]), t.φ, t.hZ, t.origin, t.rotation)
    outer_mantle = FullConeMantle{T,:inwards}( (t.r[1], t.r[2][2]), t.φ, t.hZ, t.origin, t.rotation)
    # e_bot = Annulus{T}(r = (t.r[1][1], t.r[1][2]), φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = EllipticalSurface{T}(r = (t.r[2][1], t.r[2][2]), φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    # normals of the surfaces show inside the volume primitives. 
    e_top, inner_mantle, outer_mantle
end

function surfaces(t::PartialBottomClosedTube{T,ClosedPrimitive}) where {T}
    # bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    inner_mantle = PartialConeMantle{T,:outwards}((t.r[1], t.r[2][1]), t.φ, t.hZ, t.origin, t.rotation)
    outer_mantle = PartialConeMantle{T,:inwards}( (t.r[1], t.r[2][2]), t.φ, t.hZ, t.origin, t.rotation)
    # e_bot = PartialAnnulus{T}(r = (t.r[1][1], t.r[1][2]), φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = EllipticalSurface{T}(r = (t.r[2][1], t.r[2][2]), φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    poly_l = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(t.r[2][1], zero(T), +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][2], zero(T), +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(   t.r[1], zero(T), -t.hZ)) )), t)
    poly_r = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(t.r[2][1], t.φ, +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(   t.r[1], t.φ, -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][2], t.φ, +t.hZ)) )), t)
    # normals of the surfaces show inside the volume primitives. 
    e_top, inner_mantle, outer_mantle, poly_l, poly_r
end


####################################################################
####################################################################


function Geometry(::Type{T}, t::Type{Cone}, dict::AbstractDict, input_units::NamedTuple, transformations::Transformations{T}) where {T}
    length_unit = input_units.length
    angle_unit = input_units.angle
    origin = get_origin(T, dict, length_unit)
    rotation = get_rotation(T, dict, angle_unit)

    r = parse_r_of_primitive(T, dict, length_unit, Cone)
    φ_interval = parse_φ_of_primitive(T, dict, angle_unit)
    φ = if φ_interval isa Tuple{T,T}
        rotation = rotation * RotZ{T}(φ_interval[1])
        φ_interval[2] - φ_interval[1]
    elseif isnothing(φ_interval)
        nothing
    else
        throw(ConfigFileError("Error when trying to parse φ from configuration file."))
    end

    hZ = if haskey(dict, "h")
        _parse_value(T, dict["h"], length_unit) / 2
    end

    # TODO: throw error if hZ == 0 and different radii at top and bottom ?

    cone = Cone{T}(ClosedPrimitive; 
        r = r, 
        φ = φ, 
        hZ = hZ, 
        origin = origin,
        rotation = rotation
    )
    return transform(cone, transformations)
end


function Dictionary(c::Cone{T, <:Any, TR})::OrderedDict{String, Any} where {T, TR}
    dict = OrderedDict{String, Any}()
    name  = "tube"
    if TR <: Real
        dict["r"] = c.r
    elseif TR <: Tuple{<:Real, <:Real}
        dict["r"] = OrderedDict{String, Any}("from" => c.r[1], "to" => c.r[2])
    else
        name = "cone"
        dict["r"] = OrderedDict{String, Any}()
        dict["r"]["bottom"] = if c.r[1] isa Tuple{<:Real} 
            c.r[1][1]
        elseif isnothing(c.r[1])
            0
        else
            OrderedDict{String, Any}("from" => isnothing(c.r[1][1]) ? 0 : c.r[1][1], "to" => c.r[1][2])
        end
        dict["r"]["top"] = if c.r[2] isa Tuple{<:Real} 
            c.r[2][1]
        elseif isnothing(c.r[2])
            0
        else 
            OrderedDict{String, Any}("from" => isnothing(c.r[2][1]) ? 0 : c.r[2][1], "to" => c.r[2][2])
        end
    end
    if !isnothing(c.φ) dict["phi"] = OrderedDict("from" => "0°", "to" => string(rad2deg(c.φ))*"°") end
    dict["h"] = 2*c.hZ
    if !iszero(c.origin) dict["origin"] = Dictionary(c.origin) end
    if !isone(c.rotation)
        d = Dictionary(c.rotation)
        if unique(keys(d)) == ["Z"]
            φ0 = mod2pi(_parse_value(T, d["Z"], internal_angle_unit))
            dict["phi"] = OrderedDict("from" => string(rad2deg(φ0))*"°", "to" => string(rad2deg(φ0 + c.φ))*"°")
        else
            dict["rotation"] = d
        end
    end
    OrderedDict{String, Any}(name => dict)
end

extremum(c::Cone{T}) where {T} = hypot(c.hZ, max((c.r...)...))