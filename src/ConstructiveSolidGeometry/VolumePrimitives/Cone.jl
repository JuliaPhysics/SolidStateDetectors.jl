"""
    struct Cone{T,CO,RT,TP} <: AbstractVolumePrimitive{T, CO}

T: Type of values, e.g. Float64
CO: ClosedPrimitive or OpenPrimitive <-> whether surface belongs to it or not

* `r::TR`: 
    * TR = Real -> Cylinder
    * TR = (Real, Real) -> Tube (r_in = r[1], r_out = r[2])
    * TR = ((Real,), (Real,)) Solid widening Cylinder  -> (r_bot = r[1][1], r_top = r[1][2])
    * TR = ((Real,Real), (Real,Real)) Solid widening Tube ->\n(r_bot_in = r[1][1], r_bot_out = r[1][2], r_top_in = r[2][1], r_top_out = r[2][2])
    * TR = (Nothing, (Real,Real)) Cone ->\n(r_bot_in = r_bot_out = 0, r_top_in = r[2][1], r_top_out = r[2][2])
    * TR = ((Real,Real), Nothing) Cone ->\n(r_bot_in = r[1][1], r_bot_out = r[1][2], r_top_in = r_top_out = 0)
    * ... (+ elliptical cases -> (a, b))
    * Not all are implemented yet

* `φ::TP`: 
    * TP = Nothing <-> Full in φ
    * ...
* `zH::T`: half hight/length of the cone
"""
@with_kw struct Cone{T,CO,RT,TP} <: AbstractVolumePrimitive{T, CO}
    r::RT = 1
    φ::TP = nothing
    hZ::T = 1

    origin::CartesianPoint{T} = zero(CartesianPoint{T})
    rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})
end

Cone{T,CO,RT,TP}( c::Cone{T,CO,RT,TP}; COT = CO,
            origin::CartesianPoint{T} = c.origin,
            rotation::SMatrix{3,3,T,9} = c.rotation) where {T,CO<:Union{ClosedPrimitive, OpenPrimitive},RT,TP} =
    Cone{T,COT,RT,TP}(c.r, c.φ, c.hZ, origin, rotation)

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
const PartialCylinder{T,CO} = Cone{T,CO,T,Tuple{T,T}} 

### Cylinder

function _in(pt::CartesianPoint, c::Cylinder{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    az = abs(pt.z) 
    az <= c.hZ + csgtol && begin
        r = hypot(pt.x, pt.y) 
        r <= c.r + csgtol
    end
end
function _in(pt::CartesianPoint, c::Cylinder{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) < c.hZ - csgtol &&
    hypot(pt.x, pt.y) < c.r - csgtol
end

function surfaces(t::Cylinder{T}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    mantle = ConeMantle{T,Tuple{T,T},Nothing,:inwards}((t.r, t.r), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = CircularArea{T}(r = t.r, φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = CircularArea{T}(r = t.r, φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    # normals of the surfaces show inside the volume primitives. 
    e_top, e_bot, mantle
end

### PartialCylinder

function _in(pt::CartesianPoint, c::PartialCylinder{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    az = abs(pt.z) 
    az <= c.hZ + csgtol && begin
        r = hypot(pt.x, pt.y) 
        r <= c.r + csgtol && 
        _in_angular_interval_closed(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)
    end
end
function _in(pt::CartesianPoint, c::PartialCylinder{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) + csgtol < c.hZ && begin
        r = hypot(pt.x, pt.y) 
        csgtol + r < c.r &&
        _in_angular_interval_open(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)
    end
end

function surfaces(t::PartialCylinder{T}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    mantle = PartialConeMantle{T,:inwards}((t.r, t.r), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = PartialCircularArea{T}(r = t.r, φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = PartialCircularArea{T}(r = t.r, φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    poly_l = _transform_into_global_coordinate_system(Quadrangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[1], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[1], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(    t.r,  t.φ[1], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(    t.r,  t.φ[1], -t.hZ)) )), t)
    poly_r = _transform_into_global_coordinate_system(Quadrangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[2], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(    t.r,  t.φ[2], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(    t.r,  t.φ[2], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[2], +t.hZ)) )), t)
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
const VaryingCylinder{T,CO} = Cone{T,CO,Tuple{T,T},Nothing} # Full in φ
const PartialVaryingCylinder{T,CO} = Cone{T,CO,Tuple{T,T},Tuple{T,T}}

function _in(pt::CartesianPoint, c::VaryingCylinder{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    az = abs(pt.z) 
    az <= c.hZ + csgtol && begin
        r = hypot(pt.x, pt.y) 
        rz = radius_at_z(c.hZ, c.r[1], c.r[2], pt.z)
        r <= rz + csgtol
    end
end
function _in(pt::CartesianPoint, c::VaryingCylinder{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) + csgtol < c.hZ &&
    csgtol + hypot(pt.x, pt.y) < radius_at_z(c.hZ, c.r[1], c.r[2], pt.z)
end

function surfaces(t::VaryingCylinder{T}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    mantle = ConeMantle{T,Tuple{T,T},Nothing,:inwards}((t.r[1], t.r[2]), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = CircularArea{T}(r = t.r[1], φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = CircularArea{T}(r = t.r[2], φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    # normals of the surfaces show inside the volume primitives. 
    e_top, e_bot, mantle
end

function _in(pt::CartesianPoint, c::PartialVaryingCylinder{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    az = abs(pt.z) 
    az <= c.hZ + csgtol && begin
        r = hypot(pt.x, pt.y) 
        rz = radius_at_z(c.hZ, c.r[1], c.r[2], pt.z)
        r <= rz + csgtol &&
        _in_angular_interval_closed(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)
    end
end
function _in(pt::CartesianPoint, c::PartialVaryingCylinder{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) + csgtol < c.hZ && begin
        r = hypot(pt.x, pt.y) 
        csgtol + r < radius_at_z(c.hZ, c.r[1], c.r[2], pt.z) &&
        _in_angular_interval_open(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)
    end
end

function surfaces(t::PartialVaryingCylinder{T}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    mantle = ConeMantle{T,Tuple{T,T},Tuple{T,T},:inwards}((t.r[1], t.r[2]), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = PartialCircularArea{T}(r = t.r[1], φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = PartialCircularArea{T}(r = t.r[2], φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    poly_l = _transform_into_global_coordinate_system(Quadrangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[1], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[1], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}( t.r[2],  t.φ[1], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}( t.r[1],  t.φ[1], -t.hZ)) )), t)
    poly_r = _transform_into_global_coordinate_system(Quadrangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[2], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}( t.r[1],  t.φ[2], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}( t.r[2],  t.φ[2], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[2], +t.hZ)) )), t)
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
const PartialVaryingTube{T,CO} = Cone{T,CO,Tuple{Tuple{T,T},Tuple{T,T}},Tuple{T,T}}      

#(r_bot_in = r[1][1], r_bot_out = r[1][2], r_top_in = r[2][1], r_top_out = r[2][2])

function _in(pt::CartesianPoint, c::VaryingTube{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    az = abs(pt.z)
    az <= c.hZ + csgtol && begin
        r = hypot(pt.x, pt.y)
        r_in = radius_at_z(c.hZ, c.r[1][1], c.r[2][1], pt.z)
        r_out = radius_at_z(c.hZ, c.r[1][2], c.r[2][2], pt.z)
        r_in - csgtol <= r  &&
        r <= r_out + csgtol
    end
end
function _in(pt::CartesianPoint, c::VaryingTube{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) + csgtol < c.hZ &&
    csgtol + radius_at_z(c.hZ, c.r[1][1], c.r[2][1], pt.z) < hypot(pt.x, pt.y) < radius_at_z(c.hZ, c.r[1][2], c.r[2][2], pt.z) - csgtol
end

function surfaces(t::VaryingTube{T}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    inner_mantle = ConeMantle{T,Tuple{T,T},Nothing,:outwards}((t.r[1][1], t.r[2][1]), t.φ, t.hZ, t.origin, t.rotation)
    outer_mantle = ConeMantle{T,Tuple{T,T},Nothing,:inwards}( (t.r[1][2], t.r[2][2]), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = Annulus{T}(r = (t.r[1][1], t.r[1][2]), φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = Annulus{T}(r = (t.r[2][1], t.r[2][2]), φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    # normals of the surfaces show inside the volume primitives. 
    e_top, e_bot, inner_mantle, outer_mantle
end

function _in(pt::CartesianPoint, c::PartialVaryingTube{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    az = abs(pt.z)
    az <= c.hZ + csgtol && begin
        r = hypot(pt.x, pt.y)
        r_in = radius_at_z(c.hZ, c.r[1][1], c.r[2][1], pt.z)
        r_out = radius_at_z(c.hZ, c.r[1][2], c.r[2][2], pt.z)
        r_in - csgtol <= r &&
        r <= r_out + csgtol &&
        _in_angular_interval_closed(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)
    end 
end
function _in(pt::CartesianPoint, c::PartialVaryingTube{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) + csgtol < c.hZ && begin
        r = hypot(pt.x, pt.y)
        csgtol + radius_at_z(c.hZ, c.r[1][1], c.r[2][1], pt.z) < r < radius_at_z(c.hZ, c.r[1][2], c.r[2][2], pt.z) - csgtol &&
        _in_angular_interval_open(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)
    end
end

function surfaces(t::PartialVaryingTube{T}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    inner_mantle = ConeMantle{T,Tuple{T,T},Tuple{T,T},:outwards}((t.r[1][1], t.r[2][1]), t.φ, t.hZ, t.origin, t.rotation)
    outer_mantle = ConeMantle{T,Tuple{T,T},Tuple{T,T},:inwards}( (t.r[1][2], t.r[2][2]), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = PartialAnnulus{T}(r = (t.r[1][1], t.r[1][2]), φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = PartialAnnulus{T}(r = (t.r[2][1], t.r[2][2]), φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    poly_l = _transform_into_global_coordinate_system(Quadrangle{T}((
        CartesianPoint(CylindricalPoint{T}(t.r[1][1],  t.φ[1], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][1],  t.φ[1], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][2],  t.φ[1], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[1][2],  t.φ[1], -t.hZ)) )), t)
    poly_r = _transform_into_global_coordinate_system(Quadrangle{T}((
        CartesianPoint(CylindricalPoint{T}(t.r[1][1],  t.φ[2], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[1][2],  t.φ[2], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][2],  t.φ[2], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][1],  t.φ[2], +t.hZ)) )), t)
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
const PartialUpwardCone{T,CO} = Cone{T,CO,Tuple{Tuple{Nothing,T},Nothing},Tuple{T,T}} 

function _in(pt::CartesianPoint, c::UpwardCone{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    az = abs(pt.z)
    az <= c.hZ + csgtol && begin
        r = hypot(pt.x, pt.y)
        r_out = radius_at_z(c.hZ, c.r[1][2], zero(T), pt.z)
        r <= r_out + csgtol
    end 
end
function _in(pt::CartesianPoint, c::UpwardCone{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) + csgtol < c.hZ &&
    hypot(pt.x, pt.y) + csgtol < radius_at_z(c.hZ, c.r[1][2], zero(T), pt.z)
end

function surfaces(t::UpwardCone{T}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    outer_mantle = ConeMantle{T,Tuple{T,T},Nothing,:inwards}( (t.r[1][2], zero(T)), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = CircularArea{T}(r = t.r[1][2], φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    # normals of the surfaces show inside the volume primitives. 
    e_bot, outer_mantle
end

function _in(pt::CartesianPoint, c::PartialUpwardCone{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    az = abs(pt.z)
    az <= c.hZ + csgtol && begin
        r = hypot(pt.x, pt.y)
        r_out = radius_at_z(c.hZ, c.r[1][2], zero(T), pt.z)
        r <= r_out + csgtol &&
        _in_angular_interval_closed(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)
    end 
end
function _in(pt::CartesianPoint, c::PartialUpwardCone{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) + csgtol < c.hZ && begin
        r = hypot(pt.x, pt.y)
        csgtol + r < radius_at_z(c.hZ, c.r[1][2], zero(T), pt.z) &&
        _in_angular_interval_open(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)
    end
end

function surfaces(t::PartialUpwardCone{T}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    outer_mantle = ConeMantle{T,Tuple{T,T},Tuple{T,T},:inwards}( (t.r[1][2], zero(T)), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = PartialCircularArea{T}(r = t.r[1][2], φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    poly_l = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[1], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[1], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[1][2],  t.φ[1], -t.hZ)) )), t)
    poly_r = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[2], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[1][2],  t.φ[2], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[2], +t.hZ)) )), t)
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
const PartialDownwardCone{T,CO} = Cone{T,CO,Tuple{Nothing,Tuple{Nothing,T}},Tuple{T,T}} 


function _in(pt::CartesianPoint, c::DownwardCone{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    az = abs(pt.z)
    az <= c.hZ + csgtol && begin
        r = hypot(pt.x, pt.y)
        r_out = radius_at_z(c.hZ, zero(T), c.r[2][2], pt.z)
        r <= r_out + csgtol
    end 
end
function _in(pt::CartesianPoint, c::DownwardCone{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) + csgtol < c.hZ &&
    hypot(pt.x, pt.y) + csgtol < radius_at_z(c.hZ, zero(T), c.r[2][2], pt.z)
end

function surfaces(t::DownwardCone{T}) where {T}
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    outer_mantle = ConeMantle{T,Tuple{T,T},Nothing,:inwards}( (zero(T), t.r[2][2]), t.φ, t.hZ, t.origin, t.rotation)
    e_top = CircularArea{T}(r = t.r[2][2], φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    # normals of the surfaces show inside the volume primitives. 
    e_top, outer_mantle
end

function _in(pt::CartesianPoint, c::PartialDownwardCone{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    az = abs(pt.z)
    az <= c.hZ + csgtol && begin
        r = hypot(pt.x, pt.y)
        r_out = radius_at_z(c.hZ, zero(T), c.r[2][2], pt.z)
        r <= r_out + csgtol &&
        _in_angular_interval_closed(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)
    end 
end
function _in(pt::CartesianPoint, c::PartialDownwardCone{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) + csgtol < c.hZ && begin
        r = hypot(pt.x, pt.y) 
        csgtol + r < radius_at_z(c.hZ, zero(T), c.r[2][2], pt.z) &&
        _in_angular_interval_open(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)
    end
end

function surfaces(t::PartialDownwardCone{T}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    outer_mantle = ConeMantle{T,Tuple{T,T},Tuple{T,T},:inwards}( (t.r[2][2], zero(T)), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = PartialCircularArea{T}(r = t.r[2][2], φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    poly_l = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[1], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[1], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][2],  t.φ[1], -t.hZ)) )), t)
    poly_r = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[2], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][2],  t.φ[2], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(zero(T),  t.φ[2], +t.hZ)) )), t)
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
const PartialTopClosedTube{T,CO} = Cone{T,CO,Tuple{Tuple{T,T},T},Tuple{T,T}} 


function _in(pt::CartesianPoint, c::TopClosedTube{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    az = abs(pt.z)
    az <= c.hZ + csgtol && begin
        r = hypot(pt.x, pt.y)
        r_in = radius_at_z(c.hZ, c.r[1][1], c.r[2], pt.z)
        r_out = radius_at_z(c.hZ, c.r[1][2], c.r[2], pt.z)
        r_in - csgtol <= r &&
        r <= r_out + csgtol 
    end
end
function _in(pt::CartesianPoint, c::TopClosedTube{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) + csgtol < c.hZ &&
    csgtol + radius_at_z(c.hZ, c.r[1][1], c.r[2], pt.z) < hypot(pt.x, pt.y) < radius_at_z(c.hZ, c.r[1][2], c.r[2], pt.z) - csgtol
end

function surfaces(t::TopClosedTube{T}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    # top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    inner_mantle = ConeMantle{T,Tuple{T,T},Nothing,:outwards}((t.r[1][1], t.r[2]), t.φ, t.hZ, t.origin, t.rotation)
    outer_mantle = ConeMantle{T,Tuple{T,T},Nothing,:inwards}( (t.r[1][2], t.r[2]), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = Annulus{T}(r = (t.r[1][1], t.r[1][2]), φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    # e_top = Annulus{T}(r = (t.r[2][1], t.r[2][2]), φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    # normals of the surfaces show inside the volume primitives. 
    e_bot, inner_mantle, outer_mantle
end

function _in(pt::CartesianPoint, c::PartialTopClosedTube{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    az = abs(pt.z)
    az <= c.hZ + csgtol && begin
        r = hypot(pt.x, pt.y)
        r_in = radius_at_z(c.hZ, c.r[1][1], c.r[2], pt.z)
        r_out = radius_at_z(c.hZ, c.r[1][2], c.r[2], pt.z)
        r_in - csgtol <= r &&
        r <= r_out + csgtol &&
        _in_angular_interval_closed(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)
    end 
end
function _in(pt::CartesianPoint, c::PartialTopClosedTube{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) + csgtol < c.hZ && begin
        r = hypot(pt.x, pt.y)
        csgtol + radius_at_z(c.hZ, c.r[1][1], c.r[2], pt.z) < r < radius_at_z(c.hZ, c.r[1][2], c.r[2], pt.z) - csgtol &&
        _in_angular_interval_open(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)
    end
end

function surfaces(t::PartialTopClosedTube{T}) where {T}
    bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    # top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    inner_mantle = ConeMantle{T,Tuple{T,T},Tuple{T,T},:outwards}((t.r[1][1], t.r[2]), t.φ, t.hZ, t.origin, t.rotation)
    outer_mantle = ConeMantle{T,Tuple{T,T},Tuple{T,T},:inwards}( (t.r[1][2], t.r[2]), t.φ, t.hZ, t.origin, t.rotation)
    e_bot = PartialAnnulus{T}(r = (t.r[1][1], t.r[1][2]), φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    # e_top = Annulus{T}(r = (t.r[2][1], t.r[2][2]), φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    poly_l = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(t.r[1][1],  t.φ[1], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(   t.r[2],  t.φ[1], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[1][2],  t.φ[1], -t.hZ)) )), t)
    poly_r = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(t.r[1][1],  t.φ[2], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[1][2],  t.φ[2], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(   t.r[2],  t.φ[2], +t.hZ)) )), t)
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
const PartialBottomClosedTube{T,CO} = Cone{T,CO,Tuple{T,Tuple{T,T}},Tuple{T,T}} 


function _in(pt::CartesianPoint, c::BottomClosedTube{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    az = abs(pt.z)
    az <= c.hZ + csgtol && begin
        r = hypot(pt.x, pt.y)
        r_in = radius_at_z(c.hZ, c.r[1], c.r[2][1], pt.z)
        r_out = radius_at_z(c.hZ, c.r[1], c.r[2][2], pt.z)
        r_in - csgtol <= r &&
        r <= r_out + csgtol 
    end
end
function _in(pt::CartesianPoint, c::BottomClosedTube{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) + csgtol < c.hZ &&
    csgtol + radius_at_z(c.hZ, c.r[1], c.r[2][1], pt.z) < hypot(pt.x, pt.y) < radius_at_z(c.hZ, c.r[1], c.r[2][2], pt.z) - csgtol
end

function surfaces(t::BottomClosedTube{T}) where {T}
    # bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    inner_mantle = ConeMantle{T,Tuple{T,T},Nothing,:outwards}((t.r[1], t.r[2][1]), t.φ, t.hZ, t.origin, t.rotation)
    outer_mantle = ConeMantle{T,Tuple{T,T},Nothing,:inwards}( (t.r[1], t.r[2][2]), t.φ, t.hZ, t.origin, t.rotation)
    # e_bot = Annulus{T}(r = (t.r[1][1], t.r[1][2]), φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = Annulus{T}(r = (t.r[2][1], t.r[2][2]), φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    # normals of the surfaces show inside the volume primitives. 
    e_top, inner_mantle, outer_mantle
end

function _in(pt::CartesianPoint, c::PartialBottomClosedTube{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    az = abs(pt.z)
    az <= c.hZ + csgtol && begin
        r = hypot(pt.x, pt.y)
        r_in = radius_at_z(c.hZ, c.r[1], c.r[2][1], pt.z)
        r_out = radius_at_z(c.hZ, c.r[1], c.r[2][2], pt.z)
        r_in - csgtol <= r &&
        r <= r_out + csgtol &&
        _in_angular_interval_closed(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)
    end 
end
function _in(pt::CartesianPoint, c::PartialBottomClosedTube{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} 
    abs(pt.z) + csgtol < c.hZ && begin
        r = hypot(pt.x, pt.y)
        csgtol + radius_at_z(c.hZ, c.r[1], c.r[2][1], pt.z) < r < radius_at_z(c.hZ, c.r[1], c.r[2][2], pt.z) - csgtol &&
        _in_angular_interval_open(atan(pt.y, pt.x), c.φ, csgtol = csgtol / r)
    end
end

function surfaces(t::PartialBottomClosedTube{T}) where {T}
    # bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
    top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
    inner_mantle = ConeMantle{T,Tuple{T,T},Tuple{T,T},:outwards}((t.r[1], t.r[2][1]), t.φ, t.hZ, t.origin, t.rotation)
    outer_mantle = ConeMantle{T,Tuple{T,T},Tuple{T,T},:inwards}( (t.r[1], t.r[2][2]), t.φ, t.hZ, t.origin, t.rotation)
    # e_bot = PartialAnnulus{T}(r = (t.r[1][1], t.r[1][2]), φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
    e_top = PartialAnnulus{T}(r = (t.r[2][1], t.r[2][2]), φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
    poly_l = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(t.r[2][1],  t.φ[1], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][2],  t.φ[1], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(   t.r[1],  t.φ[1], -t.hZ)) )), t)
    poly_r = _transform_into_global_coordinate_system(Triangle{T}((
        CartesianPoint(CylindricalPoint{T}(t.r[2][1],  t.φ[2], +t.hZ)),
        CartesianPoint(CylindricalPoint{T}(   t.r[1],  t.φ[2], -t.hZ)),
        CartesianPoint(CylindricalPoint{T}(t.r[2][2],  t.φ[2], +t.hZ)) )), t)
    # normals of the surfaces show inside the volume primitives. 
    e_top, inner_mantle, outer_mantle, poly_l, poly_r
end


####################################################################
####################################################################


function Geometry(::Type{T}, t::Type{Cone}, dict::AbstractDict, input_units::NamedTuple, transformations::Transformations{T}) where {T}
    length_unit = input_units.length
    angle_unit = input_units.angle
    r = parse_r_of_primitive(T, dict, length_unit, Cone)
    φ = parse_φ_of_primitive(T, dict, angle_unit)
    z = parse_height_of_primitive(T, dict, length_unit)
    hZ = typeof(z) <: Real ? z : width(z)/2
    origin = if typeof(z) <: Real
        CartesianPoint{T}(zero(T), zero(T), zero(T))
    else
        CartesianPoint{T}(zero(T), zero(T), mean(z))
    end
    cone = Cone{T, ClosedPrimitive, typeof(r), typeof(φ)}(r = r, φ = φ, hZ = hZ, origin = origin)
    return transform(cone, transformations)
end

# Cylinder


# PartialCylinder






# Tube
# function _in(pt::CartesianPoint, c::Tube{T,ClosedPrimitive}) where {T} 
#     abs(pt.z) <= c.hZ &&
#     c.r[1] <= hypot(pt.x, pt.y) <= c.r[2]
# end
# function _in(pt::CartesianPoint, c::Tube{T,OpenPrimitive}) where {T} 
#     abs(pt.z) < c.hZ &&
#     c.r[1] < hypot(pt.x, pt.y) < c.r[2]
# end

# function surfaces(t::Tube{T}) where {T}
#     bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
#     top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
#     inner_mantle = ConeMantle{T}((t.r[1], t.r[1]), t.φ, t.hZ, t.origin, t.rotation)
#     outer_mantle = ConeMantle{T}((t.r[2], t.r[2]), t.φ, t.hZ, t.origin, t.rotation)
#     e_bot = Annulus{T}(r = t.r, φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
#     e_top = Annulus{T}(r = t.r, φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
#     # normals of the surfaces show inside the volume primitives. 
#     e_top, e_bot, inner_mantle, outer_mantle
# end

# PartialTube
# function _in(pt::CartesianPoint, c::PartialTube{T,ClosedPrimitive}) where {T} 
#     abs(pt.z) <= c.hZ &&
#     c.r[1] <= hypot(pt.x, pt.y) <= c.r[2] &&
#     _in_angular_interval_closed(atan(pt.y, pt.x), c.φ)
# end
# function _in(pt::CartesianPoint, c::PartialTube{T,OpenPrimitive}) where {T} 
#     abs(pt.z) < c.hZ &&
#     c.r[1] < hypot(pt.x, pt.y) < c.r &&
#     _in_angular_interval_open(atan(pt.y, pt.x), c.φ)
# end

# function surfaces(t::PartialTube{T}) where {T}
#     bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
#     top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
#     inner_mantle = PartialConeMantle{T}((t.r[1], t.r[1]), t.φ, t.hZ, t.origin, t.rotation)
#     outer_mantle = PartialConeMantle{T}((t.r[2], t.r[2]), t.φ, t.hZ, t.origin, t.rotation)
#     e_bot = PartialAnnulus{T}(r = t.r, φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
#     e_top = PartialAnnulus{T}(r = t.r, φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
#     poly_l = _transform_into_global_coordinate_system(Quadrangle{T}((
#         CartesianPoint(CylindricalPoint{T}( t.r[1], t.φ[1], -t.hZ)),
#         CartesianPoint(CylindricalPoint{T}( t.r[2], t.φ[1], -t.hZ)),
#         CartesianPoint(CylindricalPoint{T}( t.r[2], t.φ[1], +t.hZ)),
#         CartesianPoint(CylindricalPoint{T}( t.r[1], t.φ[1], +t.hZ)) )), t)
#     poly_r = _transform_into_global_coordinate_system(Quadrangle{T}((
#         CartesianPoint(CylindricalPoint{T}( t.r[2], t.φ[2], -t.hZ)),
#         CartesianPoint(CylindricalPoint{T}( t.r[2], t.φ[2], +t.hZ)),
#         CartesianPoint(CylindricalPoint{T}( t.r[1], t.φ[2], +t.hZ)),
#         CartesianPoint(CylindricalPoint{T}( t.r[1], t.φ[2], -t.hZ)) )), t)
#     # normals of the surfaces show inside the volume primitives. 
#     e_top, e_bot, inner_mantle, outer_mantle, poly_l, poly_r
# end

# VaryingTube
# function _in(pt::CartesianPoint, c::VaryingTube{T,ClosedPrimitive}) where {T} 
#     abs(pt.z) <= c.hZ &&
#     c.r[1] <= hypot(pt.x, pt.y) <= c.r[2]
# end
# function _in(pt::CartesianPoint, c::VaryingTube{T,OpenPrimitive}) where {T} 
#     abs(pt.z) < c.hZ &&
#     c.r[1] < hypot(pt.x, pt.y) < c.r[2]
# end

# function surfaces(t::VaryingTube{T}) where {T}
#     bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
#     top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
#     inner_mantle = CylinderMantle{T}(t.r[1], t.φ, t.hZ, t.origin, t.rotation)
#     outer_mantle = CylinderMantle{T}(t.r[2], t.φ, t.hZ, t.origin, t.rotation)
#     e_bot = Annulus{T}(r = t.r, φ = t.φ, origin = bot_center_pt, rotation = t.rotation)
#     e_top = Annulus{T}(r = t.r, φ = t.φ, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
#     # normals of the surfaces show inside the volume primitives. 
#     e_top, e_bot, inner_mantle, outer_mantle
# end


# function surfaces(t::VaryingTube{T,CO}) where {T,CO}
#     bot_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -t.hZ), t) 
#     top_center_pt = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +t.hZ), t) 
#     in_mantle  = ConeMantle{T,Tuple{T,T},Nothing}((t.r[1][1], t.r[2][1]), t.φ, t.hZ, t.origin, t.rotation)
#     out_mantle = ConeMantle{T,Tuple{T,T},Nothing}((t.r[1][2], t.r[2][2]), t.φ, t.hZ, t.origin, t.rotation)
#     e_bot = EllipticalSurface{T,Tuple{T,T},Nothing}(r = t.r[1], φ = nothing, origin = bot_center_pt, rotation = t.rotation)
#     e_top = EllipticalSurface{T,Tuple{T,T},Nothing}(r = t.r[2], φ = nothing, origin = top_center_pt, rotation = -t.rotation * RotZ{T}(π))
#     # normals of the surfaces show inside the volume primitives. 
#     e_top, e_bot, in_mantle, out_mantle
# end

# PartialVaryingTube





################################
################################
################################


# function Dictionary(g::Cone{T,<:Union{T, AbstractInterval}}) where {T}
#     dict = OrderedDict{String,Any}()
#     dict["r"] = typeof(g.r) == T ? g.r : OrderedDict{String,Any}("from" => g.r.left, "to" => g.r.right)
#     if !isnothing(g.φ) dict["phi"] = OrderedDict{String,Any}("from" => g.φ.left, "to" => g.φ.right) end
#     typeof(g.z) == T ? dict["h"] = g.z * 2 : dict["z"] = OrderedDict{String,Any}("from" => g.z.left, "to" => g.z.right) 
#     OrderedDict{String,Any}("tube" => dict)
# end

# function Dictionary(g::Cone{T,<:Tuple}) where {T}
#     dict = OrderedDict{String,Any}()
#     dict["r"] = OrderedDict{String,Any}()
#     dict["r"]["bottom"] = typeof(g.r[1]) == T ? g.r[1] : OrderedDict{String,Any}("from" => g.r[1].left, "to" => g.r[1].right)
#     dict["r"]["top"] = typeof(g.r[2]) == T ? g.r[2] : OrderedDict{String,Any}("from" => g.r[2].left, "to" => g.r[2].right)
#     if !isnothing(g.φ) dict["phi"] = OrderedDict{String,Any}("from" => g.φ.left, "to" => g.φ.right) end
#     typeof(g.z) == T ? dict["h"] = g.z * 2 : dict["z"] = OrderedDict{String,Any}("from" => g.z.left, "to" => g.z.right) 
#     OrderedDict{String,Any}("cone" => dict)
# end


# get_r_limits(c::Cone{T, <:Union{T, AbstractInterval{T}}, <:Any, <:Any}) where {T} =
#     (_left_radial_interval(c.r),_right_radial_interval(c.r),_left_radial_interval(c.r),_right_radial_interval(c.r))
# get_r_limits(c::Cone{T, <:Tuple, <:Any, <:Any}) where {T} =
#     (_left_radial_interval(c.r[1]),_right_radial_interval(c.r[1]),_left_radial_interval(c.r[2]),_right_radial_interval(c.r[2]))

# get_φ_limits(c::Cone{T, <:Any, Nothing, <:Any}) where {T} = (T(0), T(2π), true)
# get_φ_limits(c::Cone{T, <:Any, <:AbstractInterval, <:Any}) where {T} = (c.φ.left, c.φ.right, false)

# get_z_limits(c::Cone{T}) where {T} = (_left_linear_interval(c.z), _right_linear_interval(c.z))

# function _is_cone_collapsed(rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T, zMin::T, zMax::T) where {T}
#     tol = geom_atol_zero(T)
#     (isapprox(rbotMin, rbotMax, atol = tol) && isapprox(rtopMin, rtopMax, atol = tol)) || isapprox(zMin, zMax, atol = tol)
# end

# function _get_decomposed_surfaces_cone(c::Cone{T}, rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T, zMin::T, zMax::T) where {T}
#     surfaces = AbstractSurfacePrimitive[]
#     #top and bottom annulus
#     tol = geom_atol_zero(T)
#     if !isapprox(rbotMin, rbotMax, atol = tol)
#         push!(surfaces, CylindricalAnnulus(c, z = zMin))
#     end
#     if !isapprox(zMin, zMax, atol = tol)
#         if !isapprox(rtopMin, rtopMax, atol = tol)
#             push!(surfaces, CylindricalAnnulus(c, z = zMax))
#         end
#         #outer conemantle
#         push!(surfaces, ConeMantle(c, rbot = rbotMax, rtop = rtopMax))
#     end
#     surfaces
# end

# #2π Cones
# function get_decomposed_surfaces(c::Cone{T, <:Union{T, Tuple{T,T}}, Nothing, <:Any}) where {T}
#     rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = get_r_limits(c)
#     zMin::T, zMax::T = get_z_limits(c)
#     _get_decomposed_surfaces_cone(c, rbotMin, rbotMax, rtopMin, rtopMax, zMin, zMax)
# end

# function get_decomposed_surfaces(c::Cone{T, <:Union{<:AbstractInterval{T}, Tuple{I,I}}, Nothing, <:Any}) where {T, I<:AbstractInterval{T}}
#     rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = get_r_limits(c)
#     zMin::T, zMax::T = get_z_limits(c)
#     surfaces = _get_decomposed_surfaces_cone(c, rbotMin, rbotMax, rtopMin, rtopMax, zMin, zMax)
#     if !_is_cone_collapsed(rbotMin, rbotMax, rtopMin, rtopMax, zMin, zMax)
#         push!(surfaces, ConeMantle(c, rbot = rbotMin, rtop = rtopMin))
#     end
#     surfaces
# end

# #non 2π Cones
# function get_decomposed_surfaces(c::Cone{T, <:Union{T, Tuple{T,T}}, <:AbstractInterval{T}, <:Any}) where {T}
#     rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = get_r_limits(c)
#     zMin::T, zMax::T = get_z_limits(c)
#     φMin::T, φMax::T, _ = get_φ_limits(c)
#     surfaces = _get_decomposed_surfaces_cone(c, rbotMin, rbotMax, rtopMin, rtopMax, zMin, zMax)
#     if !_is_cone_collapsed(rbotMin, rbotMax, rtopMin, rtopMax, zMin, zMax)
#         push!(surfaces, ConalPlane(c, φ = φMin), ConalPlane(c, φ = φMax))
#     end
#     surfaces
# end

# function get_decomposed_surfaces(c::Cone{T, <:Union{<:AbstractInterval{T}, Tuple{I,I}}, <:AbstractInterval{T}, <:Any}) where {T, I<:AbstractInterval{T}}
#     rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = get_r_limits(c)
#     zMin::T, zMax::T = get_z_limits(c)
#     φMin::T, φMax::T, _ = get_φ_limits(c)
#     surfaces = _get_decomposed_surfaces_cone(c, rbotMin, rbotMax, rtopMin, rtopMax, zMin, zMax)
#     if !_is_cone_collapsed(rbotMin, rbotMax, rtopMin, rtopMax, zMin, zMax)
#         push!(surfaces, ConalPlane(c, φ = φMin), ConalPlane(c, φ = φMax))
#         push!(surfaces, ConeMantle(c, rbot = rbotMin, rtop = rtopMin))
#     end
#     surfaces
# end

# function sample(c::Cone{T}, step::Real)::Vector{CylindricalPoint{T}} where {T}
#     zMin::T, zMax::T = get_z_limits(c)
#     φMin::T, φMax::T, _ = get_φ_limits(c)
#     samples = [
#         CylindricalPoint{T}(r,φ,z)
#         for z in zMin:step:zMax
#         for r in _left_radial_interval(get_r_at_z(c, z)):step:_right_radial_interval(get_r_at_z(c, z))
#         for φ in (r == 0 ? φMin : φMin:step/r:φMax)
#     ]
# end

# function sample(c::Cone{T}, Nsamps::NTuple{3,Int})::Vector{CylindricalPoint{T}} where {T}
#     zMin::T, zMax::T = get_z_limits(c)
#     φMin::T, φMax::T, _ = get_φ_limits(c)
#     samples = [
#         CylindricalPoint{T}(r,φ,z)
#         for z in (Nsamps[3] ≤ 1 ? zMin : range(zMin, zMax, length = Nsamps[3]))
#         for r in (Nsamps[1] ≤ 1 ? _left_radial_interval(get_r_at_z(c, z)) : range(_left_radial_interval(get_r_at_z(c, z)), _right_radial_interval(get_r_at_z(c, z)), length = Nsamps[1]))
#         for φ in (Nsamps[2] ≤ 1 ? φMin : range(φMin, φMax, length = Nsamps[2]))
#     ]
# end


#=
    This sampling might be better to be done more general:
    kinda like: sample.(lines.(surfaces.(p::VolumePrimitive))
    So we only have to right a set of sample methods for certain line types
=#
# function sample(t::Cylinder{T}; n = 4) where {T}
#     # this could be improved performance-vise, 
#     # but not that important right now as it is only called 
#     # in the initzialiaton of the grid. 
#     ehZ = CartesianPoint(zero(T), zero(T), t.hZ)
#     e_bot = Ellipse(t.r, t.φ, t.origin - t.rotation * ehZ, t.rotation)
#     e_top = Ellipse(t.r, t.φ, t.origin + t.rotation * ehZ, t.rotation)
#     φs = range(0, step = 2π / n, length = n)
#     pts_bot = [CartesianPoint(CylindricalPoint{T}(e_bot.r, φ, zero(T))) for φ in φs]
#     pts_bot = map(p -> _transform_into_global_coordinate_system(p, e_bot), pts_bot)
#     pts_top = [CartesianPoint(CylindricalPoint{T}(e_top.r, φ, zero(T))) for φ in φs]
#     pts_top = map(p -> _transform_into_global_coordinate_system(p, e_top), pts_top)
#     vcat(pts_bot, pts_top)
# end
# function sample(t::Tube{T}; n = 4) where {T}
#     ehZ = CartesianPoint(zero(T), zero(T), t.hZ)
#     φs = range(0, step = 2π / n, length = n)
#     in_e_bot  = Ellipse(t.r[1], t.φ, t.origin - t.rotation * ehZ, t.rotation)
#     out_e_bot = Ellipse(t.r[2], t.φ, t.origin - t.rotation * ehZ, t.rotation)
#     in_e_top  = Ellipse(t.r[1], t.φ, t.origin + t.rotation * ehZ, t.rotation)
#     out_e_top = Ellipse(t.r[2], t.φ, t.origin + t.rotation * ehZ, t.rotation)
#     pts_in_bot  = [CartesianPoint(CylindricalPoint{T}(in_e_bot.r, φ, zero(T))) for φ in φs]
#     pts_out_bot = [CartesianPoint(CylindricalPoint{T}(out_e_bot.r, φ, zero(T))) for φ in φs]
#     pts_in_top  = [CartesianPoint(CylindricalPoint{T}(in_e_top.r, φ, zero(T))) for φ in φs]
#     pts_out_top = [CartesianPoint(CylindricalPoint{T}(out_e_top.r, φ, zero(T))) for φ in φs]
#     pts_in_bot  = map(p -> _transform_into_global_coordinate_system(p, in_e_bot), pts_in_bot)
#     pts_out_bot = map(p -> _transform_into_global_coordinate_system(p, out_e_bot), pts_out_bot)
#     pts_in_top  = map(p -> _transform_into_global_coordinate_system(p, in_e_top), pts_in_top)
#     pts_out_top = map(p -> _transform_into_global_coordinate_system(p, out_e_top), pts_out_top)
#     vcat(pts_in_bot, pts_out_bot, pts_in_top, pts_out_top)
# end
# function sample(t::Cone{T,CO,Tuple{Tuple{T,T},Tuple{T,T}}}; n = 4) where {T,CO}
#     ehZ = CartesianPoint(zero(T), zero(T), t.hZ)
#     e_in_bot  = Ellipse(t.r[1][1], t.φ, t.origin - t.rotation * ehZ, t.rotation)
#     e_out_bot = Ellipse(t.r[1][2], t.φ, t.origin - t.rotation * ehZ, t.rotation)
#     e_in_top  = Ellipse(t.r[2][1], t.φ, t.origin + t.rotation * ehZ, t.rotation)
#     e_out_top = Ellipse(t.r[2][2], t.φ, t.origin + t.rotation * ehZ, t.rotation)
#     φs = range(0, step = 2π / n, length = n)
#     pts_in_bot = [CartesianPoint(CylindricalPoint{T}(e_in_bot.r, φ, zero(T))) for φ in φs]
#     pts_in_bot = map(p -> _transform_into_global_coordinate_system(p, e_in_bot), pts_in_bot)
#     pts_out_bot = [CartesianPoint(CylindricalPoint{T}(e_out_bot.r, φ, zero(T))) for φ in φs]
#     pts_out_bot = map(p -> _transform_into_global_coordinate_system(p, e_out_bot), pts_out_bot)
#     pts_in_top = [CartesianPoint(CylindricalPoint{T}(e_in_top.r, φ, zero(T))) for φ in φs]
#     pts_in_top = map(p -> _transform_into_global_coordinate_system(p, e_in_top), pts_in_top)
#     pts_out_top = [CartesianPoint(CylindricalPoint{T}(e_out_top.r, φ, zero(T))) for φ in φs]
#     pts_out_top = map(p -> _transform_into_global_coordinate_system(p, e_out_top), pts_out_top)
#     vcat(pts_in_bot, pts_out_bot, pts_in_top, pts_out_top)
# end

# @inline sample(c::Cone{T}) where {T} = sample(c, (2,3,3))
# @inline sample(c::Cone{T, <:Any, Nothing}) where {T} = sample(c, (2,5,3))
