"""
    struct EllipticalSurface{T,TR,TP} <: AbstractPlanarSurfacePrimitive{T}

Surface primitive describing circular bases, e.g. the top or bottom base of a [`Cone`](@ref).

## Parametric types
* `T`: Precision type.
* `TR`: Type of the radius `r`.
    * `TR == T`: Full Circle (constant radius `r`, no cut-out).
    * `TR == Tuple{T, T}`: Circular Annulus (inner radius at `r[1]`, outer radius at `r[2]`).
* `TP`: Type of the angular range `φ`.
    * `TP == Nothing`: Full 2π Cone.
    * `TP == Tuple{T, T}`: Partial Cone ranging from `φ[1]` to `φ[2]`.
    
## Fields
* `r::TR`: Definition of the radius of the `EllipticalSurface` (in m).
* `φ::TP`: Range in polar angle `φ` over which the `EllipticalSurface` extends (in radians).
* `origin::CartesianPoint{T}`: The position of the center of the `EllipticalSurface`.
* `rotation::SMatrix{3,3,T,9}`: Matrix that describes a rotation of the `EllipticalSurface` around its `origin`.
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

normal(es::EllipticalSurface{T}, ::CartesianPoint{T} = zero(CartesianPoint{T})) where {T} = es.rotation * CartesianVector{T}(zero(T), zero(T), one(T))

extremum(es::EllipticalSurface{T,T}) where {T} = es.r
extremum(es::EllipticalSurface{T,Tuple{T,T}}) where {T} = es.r[2] # r_out always larger r_in: es.r[2] > es.r[2]

function lines(sp::CircularArea{T}; n = 2) where {T} 
    circ = Circle{T}(r = sp.r, φ = sp.φ, origin = sp.origin, rotation = sp.rotation)
    φs = range(T(0), step = T(2π) / n, length = n)
    edges = [ Edge{T}(_transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(zero(T), φ, zero(T))), sp),
                      _transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r, φ, zero(T))), sp)) for φ in φs ]
    return (circ, edges)
end
function lines(sp::PartialCircularArea{T}; n = 2) where {T} 
    circ = PartialCircle{T}(r = sp.r, φ = sp.φ, origin = sp.origin, rotation = sp.rotation)
    φs = range(sp.φ[1], stop = sp.φ[2], length = n)
    edges = [ Edge{T}(_transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(zero(T), φ, zero(T))), sp),
                      _transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r, φ, zero(T))), sp)) for φ in φs ]
    return (circ, edges)
end

function lines(sp::Annulus{T}; n = 2) where {T} 
    circ_in  = Circle{T}(r = sp.r[1], φ = sp.φ, origin = sp.origin, rotation = sp.rotation)
    circ_out = Circle{T}(r = sp.r[2], φ = sp.φ, origin = sp.origin, rotation = sp.rotation)
    φs = range(T(0), step = T(2π) / n, length = n)
    edges = [ Edge{T}(_transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[1], φ, zero(T))), sp),
                      _transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[2], φ, zero(T))), sp)) for φ in φs ]
    return (circ_in, circ_out, edges)
end
function lines(sp::PartialAnnulus{T}; n = 2) where {T} 
    circ_in  = PartialCircle{T}(r = sp.r[1], φ = sp.φ, origin = sp.origin, rotation = sp.rotation)
    circ_out = PartialCircle{T}(r = sp.r[2], φ = sp.φ, origin = sp.origin, rotation = sp.rotation)
    φs = range(sp.φ[1], stop = sp.φ[2], length = n)
    edges = [ Edge{T}(_transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[1], φ, zero(T))), sp),
                      _transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[2], φ, zero(T))), sp)) for φ in φs ]
    return (circ_in, circ_out, edges)
end

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
