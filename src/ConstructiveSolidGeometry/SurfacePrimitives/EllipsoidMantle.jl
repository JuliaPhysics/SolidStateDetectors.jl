@with_kw struct EllipsoidMantle{T,RT,PT,TT} <: AbstractCurvedSurfacePrimitive{T}
    r::RT = 1
    φ::PT = nothing
    θ::TT = nothing

    origin::CartesianPoint{T} = zero(CartesianPoint{T})
    rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})
end

const FullSphereMantle{T} = EllipsoidMantle{T,T,Nothing,Nothing}
const FullEllipsoidMantle{T} = EllipsoidMantle{T,NTuple{3,T},Nothing,Nothing}

function lines(em::FullSphereMantle{T}) where {T} 
    ellipse_xy = Ellipse{T,T,Nothing}(r = em.r[1], φ = em.φ, origin = em.origin, rotation = em.rotation)
    ellipse_xz = Ellipse{T,T,Nothing}(r = em.r[1], φ = em.φ, origin = em.origin, rotation = em.rotation * RotX(T(π)/2))
    ellipse_yz = Ellipse{T,T,Nothing}(r = em.r[1], φ = em.φ, origin = em.origin, rotation = em.rotation * RotX(T(π)/2) * RotY(T(π)/2))
    (ellipse_xy, ellipse_xz, ellipse_yz)
end
function lines(em::FullEllipsoidMantle{T}) where {T} 
    ellipse_xy = Ellipse{T,NTuple{2, Tuple{T}},Nothing}(r = ((em.r[1],), (em.r[2],)), φ = em.φ, origin = em.origin, rotation = em.rotation)
    ellipse_xz = Ellipse{T,NTuple{2, Tuple{T}},Nothing}(r = ((em.r[1],), (em.r[3],)), φ = em.φ, origin = em.origin, rotation = em.rotation * RotX(T(π)/2))
    ellipse_yz = Ellipse{T,NTuple{2, Tuple{T}},Nothing}(r = ((em.r[2],), (em.r[3],)), φ = em.φ, origin = em.origin, rotation = em.rotation * RotX(T(π)/2) * RotY(T(π)/2))
    (ellipse_xy, ellipse_xz, ellipse_yz)
end

extremum(e::EllipsoidMantle{T,T}) where {T} = e.r
extremum(e::EllipsoidMantle{T,NTuple{3,T}}) where {T} = max(e.r)

"""
    intersection(cm::EllipsoidMantle{T,NTuple{3,T}}, l::Line{T}) where {T}

The function will always return 2 CartesianPoint's.
If the line just touches the mantle, the two points will be the same. 
If the line does not touch the mantle at all, the two points will have NaN's as there coordinates.
"""
function intersection(em::EllipsoidMantle{T,NTuple{3,T}}, l::Line{T}) where {T}
    obj_l = _transform_into_object_coordinate_system(l, em) # direction is not normalized
    
    L1 = obj_l.origin.x
    L2 = obj_l.origin.y
    L3 = obj_l.origin.z
    D1 = obj_l.direction.x
    D2 = obj_l.direction.y
    D3 = obj_l.direction.z
    
    R1 = em.r[1]
    R2 = em.r[2]
    R3 = em.r[3]

    term3 = D1^2/R1^2 + D2^2/R2^2 + D3^2/R3^2
    term1 = ((2D1*L1)/R1^2 + (2D2*L2)/R2^2 + (2D3*L3)/R3^2)^2 - 4*term3*(L1^2/R1^2 + L2^2/R2^2 + L3^2/R3^2 - 1)
    term2 = -(2D1*L1)/R1^2 - (2D2*L2)/R2^2 - (2D3*L3)/R3^2

    term1 = term1 < 0 ? sqrt(abs(term1)) : sqrt(term1)

    λ1 = (term2 - term1) / 2term3
    λ2 = (term2 + term1) / 2term3

    ints1 = obj_l.origin + λ1 * obj_l.direction 
    ints2 = obj_l.origin + λ2 * obj_l.direction 
    return _transform_into_global_coordinate_system(ints1, em), 
           _transform_into_global_coordinate_system(ints2, em)
end
"""
    intersection(cm::EllipsoidMantle{T,T}, l::Line{T}) where {T}

The function will always return 2 CartesianPoint's.
If the line just touches the mantle, the two points will be the same. 
If the line does not touch the mantle at all, the two points will have NaN's as there coordinates.
"""
function intersection(em::EllipsoidMantle{T,T}, l::Line{T}) where {T}
    obj_l = _transform_into_object_coordinate_system(l, em) # direction is not normalized
    
    L1 = obj_l.origin.x
    L2 = obj_l.origin.y
    L3 = obj_l.origin.z
    D1 = obj_l.direction.x
    D2 = obj_l.direction.y
    D3 = obj_l.direction.z
    
    R = em.r

    term3 = D1^2 + D2^2 + D3^2
    term1 = (2D1*L1 + 2D2*L2 + 2D3*L3)^2 - 4*term3*(L1^2 + L2^2 + L3^2 - R^2)
    term2 = -D1*L1 - D2*L2 - D3*L3
    
    if term1 < 0 term1 = abs(term1) end

    λ1 = (term2 - sqrt(term1)/2) / term3
    λ2 = (term2 + sqrt(term1)/2) / term3

    ints1 = obj_l.origin + λ1 * obj_l.direction 
    ints2 = obj_l.origin + λ2 * obj_l.direction 
    return _transform_into_global_coordinate_system(ints1, em), 
           _transform_into_global_coordinate_system(ints2, em)
end

# x = (-1/2 sqrt((2 D1 L1 + 2 D2 L2 + 2 D3 L3)^2 - 4 (D1^2 + D2^2 + D3^2) (L1^2 + L2^2 + L3^2 - R^2)) 
#     - D1 L1 - D2 L2 - D3 L3)/(D1^2 + D2^2 + D3^2) 
#     #and D1^2 + D2^2 + D3^2!=0 and R!=0

# y = (+1/2 sqrt((2 D1 L1 + 2 D2 L2 + 2 D3 L3)^2 - 4 (D1^2 + D2^2 + D3^2) (L1^2 + L2^2 + L3^2 - R^2)) 
#     - D1 L1 - D2 L2 - D3 L3)/(D1^2 + D2^2 + D3^2) 
    #and D1^2 + D2^2 + D3^2!=0 and R!=0



# get_r_limits(s::SphereMantle{T}) where {T} = (_left_radial_interval(s.r), _right_radial_interval(s.r))
# get_φ_limits(s::SphereMantle{T}) where {T} = (T(0), T(2π), true)
# get_z_limits(s::SphereMantle{T}) where {T} = (-_right_linear_interval(s.r), _right_linear_interval(s.r))

# in(p::AbstractCoordinatePoint, s::SphereMantle) = _isapprox_sph_r(p, s.r)

# #=
# function sample(s::SphereMantle{T}, step::Real)::Vector{CylindricalPoint{T}} where {T}
#     samples = [
#         CylindricalPoint{T}(sqrt(s.r^2 - z^2),φ,z)
#         for z in in -s.r:step:s.r
#         for φ in 0:step/s.r:2π
#     ]
# end
# =#

# function sample(s::SphereMantle{T}, Nsamps::NTuple{3,Int})::Vector{CylindricalPoint{T}} where {T}
#     samples = [
#         CylindricalPoint{T}(sqrt(s.r^2 - z^2),φ,z)
#         for z in (Nsamps[3] ≤ 1 ? s.r : range(-s.r, s.r, length = Nsamps[3]))
#         for φ in (Nsamps[2] ≤ 1 ? 0 : range(0, 2π, length = Nsamps[2]))
#     ]
# end

# function sample(s::SphereMantle{T}, g::CylindricalTicksTuple{T})::Vector{CylindricalPoint{T}} where {T}
#     samples = [
#         CylindricalPoint{T}((s.r^2 - z^2 < 0 ? 0 : sqrt(s.r^2 - z^2)),φ,z)
#         for z in get_z_ticks(s, g)
#         for φ in get_φ_ticks(s, g)
#     ]
# end

# function _get_x_at_z(s::SphereMantle{T}, g::CartesianTicksTuple{T}, z::T) where {T}
#     if s.r < abs(z) return [0] end
#     R::T = sqrt(s.r^2 - z^2)
#     x_from_y::Vector{T} = sqrt.(R^2 .- filter(y -> abs(y) <= R, g.y).^2)
#     _get_ticks(sort!(vcat(g.x, x_from_y, -x_from_y)), -R, R)
# end

# function _get_y_at_z(s::SphereMantle{T}, z::T, x::T) where {T}
#     tmp::T = s.r^2 - hypot(x,z)^2
#     if tmp < 0 return (0,) end
#     (-sqrt(tmp), sqrt(tmp))
# end

# function sample(s::SphereMantle{T}, g::CartesianTicksTuple{T})::Vector{CartesianPoint{T}} where {T}
#     samples = [
#         CartesianPoint{T}(x,y,z)
#         for z in get_z_ticks(s, g)
#         for x in _get_x_at_z(s, g, z)
#         for y in _get_y_at_z(s, z, x)
#     ]
# end

# distance_to_surface(point::CylindricalPoint{T}, s::SphereMantle{T}) where {T} = abs(hypot(point.r, point.z) - s.r)

# distance_to_surface(point::CartesianPoint{T}, s::SphereMantle{T}) where {T} = abs(norm(point) - s.r)
