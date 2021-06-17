struct Plane{T} <: AbstractFlatSurfacePrimitive{T}
    origin::CartesianPoint{T}
    normal::CartesianVector{T}
    Plane{T}(o, n) where {T} = new{T}(o, normalize(n))
end
Plane(origin::CartesianPoint{T}, normal::CartesianVector{T}) where {T} = Plane{T}(origin, normal)

normal(p::Plane) = p.normal
origin(p::Plane) = p.origin

isinfront(pt::AbstractCoordinatePoint, p::Plane) = (pt - origin(p)) ⋅ normal(p) > 0
isbehind(pt::AbstractCoordinatePoint, p::Plane) = (pt - origin(p)) ⋅ normal(p) < 0
in(pt::AbstractCoordinatePoint, p::Plane) = (pt - origin(p)) ⋅ normal(p) == 0

_distance(pt::AbstractCoordinatePoint, p::Plane) = (pt - origin(p)) ⋅ normal(p)
distance(pt::AbstractCoordinatePoint, p::Plane) = abs(_distance(pt, p))


function evaluate(p::Plane{T}, x, y, ::Val{:xy}) where {T}
    # plane.normal may not be perpendicular to the z-axis
    line = Line{T}(CartesianPoint{T}(x, y, zero(T)), CartesianVector{T}(zero(T), zero(T), one(T)))
    z = (p.normal ⋅ p.origin - p.normal ⋅ line.origin) / (p.normal ⋅ normalize(line.direction))
    CartesianPoint{T}(x, y, z)
end
function evaluate(p::Plane{T}, x, z, ::Val{:xz}) where {T}
    # plane.normal may not be perpendicular to the y-axis
    line = Line{T}(CartesianPoint{T}(x, zero(T), z), CartesianVector{T}(zero(T), one(T), zero(T)))
    y = (p.normal ⋅ p.origin - p.normal ⋅ line.origin) / (p.normal ⋅ normalize(line.direction))
    CartesianPoint{T}(x, y, z)
end
function evaluate(p::Plane{T}, y, z, ::Val{:yz}) where {T}
    # plane.normal may not be perpendicular to the x-axis
    line = Line{T}(CartesianPoint{T}(zero(T), y, z), CartesianVector{T}(one(T), zero(T), zero(T)))
    x = (p.normal ⋅ p.origin - p.normal ⋅ line.origin) / (p.normal ⋅ normalize(line.direction))
    CartesianPoint{T}(x, y, z)
end


function evaluate(p::Plane{T}, r, φ, ::Val{:rφ}) where {T}
    # plane.normal may not be perpendicular to the z-axis
    car_pt = CartesianPoint(CylindricalPoint{T}(r, φ, zero(T)))
    line = Line{T}(CartesianPoint{T}(car_pt.x, car_pt.y, zero(T)), CartesianVector{T}(zero(T), zero(T), one(T)))
    z = (p.normal ⋅ p.origin - p.normal ⋅ line.origin) / (p.normal ⋅ normalize(line.direction))
    CylindricalPoint(CartesianPoint{T}(car_pt.x, car_pt.y, z))
end
# function evaluate(p::Plane{T}, r, z, ::Val{:rz}) where {T}
#     We don't need this
# end
function evaluate(p::Plane{T}, φ, z, ::Val{:φz}) where {T}
    zero(CylindricalPoint{T})
end