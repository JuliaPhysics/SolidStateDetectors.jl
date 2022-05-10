"""
    struct Plane{T} <: AbstractPlanarSurfacePrimitive{T}

Surface primitive describing a two-dimensional flat plane in three-dimensional space.
    
## Fields
* `origin::CartesianPoint{T}`: Point in the `Plane`.
* `normal::CartesianVector{T}`: Normal vector of the `Plane`, normalized to length 1.
"""
struct Plane{T} <: AbstractPlanarSurfacePrimitive{T}
    origin::CartesianPoint{T}
    normal::CartesianVector{T}
    Plane{T}(o, n) where {T} = new{T}(o, normalize(n))
end

#Type promotion happens here
function Plane(origin::PT, normal::NO) where {PT,NO}
    eltypes = _csg_get_promoted_eltype.((PT,NO))
    T = float(promote_type(eltypes...))
    Plane{T}(origin, normal)
end

function Plane(;
    origin = zero(CartesianPoint{Int}), 
    normal = CartesianVector{Int}(0,0,1)
)
    Plane(origin, normal)
end

function Plane{T}(;
    origin = zero(CartesianPoint{Float64}), 
    normal = CartesianVector{Float64}(0,0,1)
) where {T}
    Plane{T}(origin, normal)
end

normal(p::Plane) = p.normal
origin(p::Plane) = p.origin

isinfront(pt::AbstractCoordinatePoint, p::Plane) = (pt - origin(p)) ⋅ normal(p) > 0
isbehind(pt::AbstractCoordinatePoint, p::Plane) = (pt - origin(p)) ⋅ normal(p) < 0

_distance(pt::AbstractCoordinatePoint, p::Plane) = (pt - origin(p)) ⋅ normal(p)
distance(pt::AbstractCoordinatePoint, p::Plane) = abs(_distance(pt, p))

"""
    intersection(p::Plane{T}, line::Line{T}) where {T}

Calculates the intersections of a `Line` with a `Plane`.

## Arguments
* `cm::Plane{T}`: The `Plane`.
* `l::Line{T}`: The `Line`.

!!! note 
    The function will always return one Point as a Tuple.
    If the line is parallel to the plane, the point will have `NaN`'s/`Inf`'s as values.
"""
function intersection(p::Plane{T}, line::Line{T}) where {T}
    ndir = normalize(line.direction)
    λ = (p.normal ⋅ p.origin - p.normal ⋅ line.origin) / (p.normal ⋅ ndir)
    (line.origin + λ * ndir,)
end

intersection(p::AbstractPlanarSurfacePrimitive{T}, line::Line{T}) where {T} = intersection(Plane(p), line)
