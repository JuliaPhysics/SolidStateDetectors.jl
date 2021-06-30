"""
    EllipticalSurface{T,TR,TP} <: AbstractSurfacePrimitive{T}

* `r::TR`: 
    * TR = Real -> Full Circle (a = b = r)
    * TR = (Real, Real) -> Circular Annulus (r_in = r[1], r_out = r[2])
    * TR = ((Real,), (Real,)) -> Full Ellipse (a = r[1][1], b = r[2][1])
    * TR = ((Real, Real),(Real, Real)) -> Elliptical Annulus \n(a_in = r[1][1], a_out = r[1][2], b_in = r[2][1], b_out = r[2][2])
    * Not all are implemented yet

* `φ::TP`: 
    * TP = Nothing <-> Full in φ
    * ...
"""
@with_kw struct EllipticalSurface{T,TR,TP} <: AbstractPlanarSurfacePrimitive{T}
    r::TR = 1
    φ::TP = nothing

    origin::CartesianPoint{T} = zero(CartesianPoint{T})
    rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})
end

const CircularArea{T} = EllipticalSurface{T,T,Nothing}
const PartialCircularArea{T} = EllipticalSurface{T,T,Tuple{T,T}}

const Annulus{T} = EllipticalSurface{T,Tuple{T,T},Nothing}
const PartialAnnulus{T} = EllipticalSurface{T,Tuple{T,T},Tuple{T,T}}

Plane(es::EllipticalSurface{T}) where {T} = Plane{T}(es.origin, es.rotation * CartesianVector{T}(zero(T),zero(T),one(T)))

normal(es::EllipticalSurface{T}) where {T} = es.rotation * CartesianVector{T}(zero(T), zero(T), one(T))

extremum(es::EllipticalSurface{T,T}) where {T} = es.r
extremum(es::EllipticalSurface{T,Tuple{T,T}}) where {T} = es.r[2] # r_out always larger r_in: es.r[2] > es.r[2]

function lines(sp::CircularArea{T}) where {T} 
    circ = Circle{T}(r = sp.r, φ = sp.φ, origin = sp.origin, rotation = sp.rotation)
    return (circ,)
end
function lines(sp::PartialCircularArea{T}) where {T} 
    circ = PartialCircle{T}(r = sp.r, φ = sp.φ, origin = sp.origin, rotation = sp.rotation)
    p_l = CartesianPoint(CylindricalPoint{T}(sp.r, sp.φ[1], zero(T)))
    p_r = CartesianPoint(CylindricalPoint{T}(sp.r, sp.φ[2], zero(T)))
    edge_l = Edge(sp.origin, _transform_into_global_coordinate_system(p_l, sp))
    edge_r = Edge(sp.origin, _transform_into_global_coordinate_system(p_r, sp))
    return (circ, edge_l, edge_r)
end

function lines(sp::Annulus{T}) where {T} 
    circ_in  = Circle{T}(r = sp.r[1], φ = sp.φ, origin = sp.origin, rotation = sp.rotation)
    circ_out = Circle{T}(r = sp.r[2], φ = sp.φ, origin = sp.origin, rotation = sp.rotation)
    return (circ_in, circ_out)
end
function lines(sp::PartialAnnulus{T}) where {T} 
    circ_in  = PartialCircle{T}(r = sp.r[1], φ = sp.φ, origin = sp.origin, rotation = sp.rotation)
    circ_out = PartialCircle{T}(r = sp.r[2], φ = sp.φ, origin = sp.origin, rotation = sp.rotation)
    p_l_in  = CartesianPoint(CylindricalPoint{T}(sp.r[1], sp.φ[1], zero(T)))
    p_l_out = CartesianPoint(CylindricalPoint{T}(sp.r[2], sp.φ[1], zero(T)))
    p_r_in  = CartesianPoint(CylindricalPoint{T}(sp.r[1], sp.φ[2], zero(T)))
    p_r_out = CartesianPoint(CylindricalPoint{T}(sp.r[2], sp.φ[2], zero(T)))
    edge_l = Edge(_transform_into_global_coordinate_system(p_l_in,  sp),
                  _transform_into_global_coordinate_system(p_l_out, sp))
    edge_r = Edge(_transform_into_global_coordinate_system(p_r_in,  sp),
                  _transform_into_global_coordinate_system(p_r_out, sp))
    return (circ_in, circ_out, edge_l, edge_r)
end

# #Constructors
# CylindricalAnnulus(c::Cone{T}; z = 0) where {T} = CylindricalAnnulus(T, get_r_at_z(c,z), c.φ, T(z))

# function CylindricalAnnulus(t::Torus{T}; θ = 0) where {T}
#     r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
#     θ = T(mod(θ,2π))
#     if θ == T(0)
#         rMin = t.r_torus + r_tubeMin
#         rMax = t.r_torus + r_tubeMax
#     elseif θ == T(π)
#         rMin = t.r_torus - r_tubeMax
#         rMax = t.r_torus - r_tubeMin
#     else
#         @error "CylindricalAnnulus not defined for torroidal cordinate θ ≠ 0 and θ ≠ π. Use ConeMantle"
#     end
#     r = rMin == 0 ? T(rMax) : T(rMin)..T(rMax)
#     CylindricalAnnulus( T, r, t.φ, t.z)
# end


# function CylindricalAnnulus(; rMin = 0, rMax = 1, φMin = 0, φMax = 2π, z = 0)
#     T = float(promote_type(typeof.((rMin, rMax, φMin, φMax, z))...))
#     r = rMin == 0 ? T(rMax) : T(rMin)..T(rMax)
#     φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
#     CylindricalAnnulus(T, r, φ, T(z))
# end

# CylindricalAnnulus(rMin, rMax, φMin, φMax, z) = CylindricalAnnulus(;rMin = rMin, rMax = rMax, φMin = φMin, φMax = φMax, z = z)

# function CylindricalAnnulus(r::Real, z::Real)
#     T = float(promote_type(typeof.((r, z))...))
#     CylindricalAnnulus(T, T(r), nothing, T(z))
# end

# get_r_limits(a::CylindricalAnnulus{T, <:Union{T, AbstractInterval{T}}, <:Any}) where {T} =
#     (_left_radial_interval(a.r), _right_radial_interval(a.r))

# get_φ_limits(a::CylindricalAnnulus{T, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
# get_φ_limits(a::CylindricalAnnulus{T, <:Any, <:AbstractInterval}) where {T} = (a.φ.left, a.φ.right, false)

# in(p::AbstractCoordinatePoint, a::CylindricalAnnulus{T, <:Any, Nothing}) where {T} = _isapprox_z(p, a.z) && _in_cyl_r(p, a.r)

# in(p::AbstractCoordinatePoint, a::CylindricalAnnulus{T, <:Any, <:AbstractInterval}) where {T} = _isapprox_z(p, a.z) && _in_φ(p, a.φ) && _in_cyl_r(p, a.r)

# #=
# function sample(a::CylindricalAnnulus{T}, step::Real)::Vector{CylindricalPoint{T}} where {T}
#     rMin::T, rMax::T = get_r_limits(a)
#     φMin::T, φMax::T, _ = get_φ_limits(a)
#     samples = [
#         CylindricalPoint{T}(r,φ,a.z)
#         for r in rMin:step:rMax
#         for φ in (r == 0 ? φMin : φMin:step/r:φMax)
#     ]
# end
# =#

# function sample(a::CylindricalAnnulus{T}, Nsamps::NTuple{3,Int})::Vector{CylindricalPoint{T}} where {T}
#     rMin::T, rMax::T = get_r_limits(a)
#     φMin::T, φMax::T, _ = get_φ_limits(a)
#     samples = [
#         CylindricalPoint{T}(r,φ,a.z)
#         for r in (Nsamps[1] ≤ 1 ? rMin : range(rMin, rMax, length = Nsamps[1]))
#         for φ in (Nsamps[2] ≤ 1 ? φMin : range(φMin, φMax, length = Nsamps[2]))
#     ]
# end

# function sample(a::CylindricalAnnulus{T}, g::CylindricalTicksTuple{T})::Vector{CylindricalPoint{T}} where {T}
#     samples = [
#         CylindricalPoint{T}(r,φ,a.z)
#         for r in get_r_ticks(a, g)
#         for φ in get_φ_ticks(a, g)
#     ]
# end

# function sample(a::CylindricalAnnulus{T}, g::CartesianTicksTuple{T})::Vector{CartesianPoint{T}} where {T}
#     L::T = _left_radial_interval(a.r)
#     R::T = _right_radial_interval(a.r)
#     samples = [
#         CartesianPoint{T}(x,y,a.z)
#         for x in _get_ticks(g.x, -R, R)
#         for y in (abs(x) > L ? _get_ticks(g.y, -sqrt(R^2 - x^2), sqrt(R^2 - x^2)) : vcat(_get_ticks(g.y, -sqrt(R^2 - x^2), -sqrt(L^2 - x^2)), _get_ticks(g.y, sqrt(L^2 - x^2), sqrt(R^2 - x^2))))
#         if a.φ === nothing || mod(atan(y, x), T(2π)) in a.φ
#     ]
# end

# function distance_to_surface(point::AbstractCoordinatePoint{T}, a::CylindricalAnnulus{T, <:Any, Nothing})::T where {T}
#     point = CylindricalPoint(point)
#     rMin::T, rMax::T = get_r_limits(a)
#     _in_cyl_r(point, a.r) ? abs(point.z - a.z) : hypot(point.z - a.z, min(abs(point.r - rMin), abs(point.r - rMax)))
# end

# function distance_to_surface(point::AbstractCoordinatePoint{T}, a::CylindricalAnnulus{T, <:Any, <:AbstractInterval})::T where {T}
#     pcy = CylindricalPoint(point)
#     rMin::T, rMax::T = get_r_limits(a)
#     φMin::T, φMax::T, _ = get_φ_limits(a)
#     if _in_φ(pcy, a.φ)
#         Δz = abs(pcy.z - a.z)
#         return _in_cyl_r(pcy, a.r) ? Δz : hypot(Δz, min(abs(pcy.r - rMin), abs(pcy.r - rMax)))
#     else
#         φNear = _φNear(pcy.φ, φMin, φMax)
#         if rMin == rMax
#             return norm(CartesianPoint(point)-CartesianPoint(CylindricalPoint{T}(rMin,φNear,a.z)))
#         else
#             return distance_to_line(CartesianPoint(point),
#                                     LineSegment(T,CartesianPoint(CylindricalPoint{T}(rMin,φNear,a.z)),
#                                                 CartesianPoint(CylindricalPoint{T}(rMax,φNear,a.z))
#                                                 )
#                                     )
#         end
#     end
# end
