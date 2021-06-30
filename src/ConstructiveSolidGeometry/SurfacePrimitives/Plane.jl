struct Plane{T} <: AbstractPlanarSurfacePrimitive{T}
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


function intersection(p::Plane{T}, line::Line{T}) where {T}
    ndir = normalize(line.direction)
    λ = (p.normal ⋅ p.origin - p.normal ⋅ line.origin) / (p.normal ⋅ ndir)
    line.origin + λ * ndir
end
