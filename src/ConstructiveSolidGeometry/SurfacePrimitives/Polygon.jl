"""
    struct Polygon{N,T} <: AbstractPlanarSurfacePrimitive{T}
        
Surface primitive describing a planar polygon, e.g. the base of a [`RegularPrism`](@ref).

## Parametric types
* `T`: Precision type.
* `N`: Number of vertices of the `Polygon`.

## Fields
* `points::SVector{N, CartesianPoint{T}}`: Vertices of the polygon in the order with which they are connected.
"""
struct Polygon{N,T} <: AbstractPlanarSurfacePrimitive{T}
    points::SVector{N, CartesianPoint{T}} 
end
const Triangle{T} = Polygon{3, T}
const Quadrangle{T} = Polygon{4, T}

normal(p::Polygon) = normalize((p.points[2] - p.points[1]) × (p.points[3] - p.points[1]))

vertices(p::Polygon) = p.points
vertices(p::Polygon, n::Int) = vertices(p)

function _sample_excluding_border(t::Triangle{T}, spacing::T)::Vector{CartesianPoint{T}} where {T}
    push = spacing/2.5#best value for not showing triangle edges
    u = t.points[2] - t.points[1]
    v = t.points[3] - t.points[1]
    w = t.points[3] - t.points[2]
    nu = norm(u)
    nv = norm(v)
    nw = norm(w)
    sθ = norm(u × v)/(nu*nv)
    sφ = norm(u × w)/(nu*nw)
    p1 = t.points[1] + (push/sθ)*(u/nu + v/nv)
    p2 = t.points[2] + (push/sφ)*(-u/nu + w/nw)
    if (p2-p1) ⋅ u > 0 
        su = norm(p2-p1)/nu
        [p1 + (su*a*u + su*b*v)
            for a in range(0, stop = 1,     length = max(2, 1 + Int(ceil(su*nu*sθ/spacing)))) 
            for b in range(0, stop = 1 - a, length = max(2, 1 + Int(ceil(su*nv*(1 - a)/spacing))))]
    else
        []
    end
end

triangles(p::Polygon{N,T}) where {N,T} = [Triangle{T}([p.points[1], p.points[i], p.points[i+1]]) for i in 2:N-1]
triangles(p::Polygon, n_arc::Int) = triangles(p)

function sample(p::Polygon{N,T}, spacing::T)::Vector{CartesianPoint{T}} where {N,T}
    v = [s for e in edges(p) for s in sample(e,n=max(2,Int(ceil(norm(e.b-e.a)/spacing))))]
    for t in triangles(p)
        append!(v, _sample_excluding_border(t, spacing))
    end
    v
end

function extremum(p::Polygon{N,T})::T where {N,T}
    c = barycenter(p.points)
    m = maximum([norm(point - c) for point in p.points])
end

connections(p::Polygon{N}) where {N} = [collect(1:N)]
connections(p::Polygon, ::Int) = connections(p)
connections(p::Polygon{N}, ::Int, ::Int) where {N} = [[i%N+1, (i+1)%N+1] for i in 0:N-1]
get_label_name(::Polygon) = "Polygon"

extreme_points(p::Polygon) = p.points

function edges(p::Polygon{N,T})::NTuple{N, Edge{T}} where {N,T}
    vs::SVector{N, CartesianPoint{T}} = vertices(p)
    NTuple{N, Edge{T}}( Edge(vs[i],vs[i%N+1]) for i in 1:N )
end

flip(p::Polygon) = Polygon(reverse(p.points))

lines(p::Polygon) = edges(p)

Plane(p::Polygon{N, T}) where {N, T} = Plane{T}(p.points[1], (p.points[2] - p.points[1]) × (p.points[3] - p.points[1]))

function _get_rot_for_rotation_on_xy_plane(p::Polygon{<:Any, T}) where {T}
    n = normal(p)
    rot_angle_rad = CartesianVector{T}(0,0,1) ⋅ n
    return if abs(rot_angle_rad) == 1
        RotationVec(zero(T), zero(T), zero(T))
    else
        rot_angle = -acos(rot_angle_rad)
        rot_axis = rot_angle * normalize(CartesianVector{T}(0,0,1) × n)
        RotationVec(rot_axis[1], rot_axis[2], rot_axis[3])
    end
end

function _rotate_on_xy_plane(p::Polygon{<:Any, T}) where {T}
    rot = _get_rot_for_rotation_on_xy_plane(p)
    map(pt -> rot * pt, p.points)
end

function in(pt::CartesianPoint{T}, p::Polygon{N,T}, csgtol::T = csg_default_tol(T)) where {N,T}
    plane = Plane(p)
    return isapprox((pt - origin(plane)) ⋅ normal(plane), 0, atol = csgtol) && begin
        rot = _get_rot_for_rotation_on_xy_plane(p)
        vs = vertices(p)

        pts2d = SVector{N+1, SVector{2,T}}(
            [SVector{2,T}(view(rot * (vs[i % N + 1] - cartesian_zero), 1:2)) for i in 0:N]...
        )

        rotated_pt = rot * (pt - cartesian_zero)
        point2d = SVector{2,T}(view(rotated_pt, 1:2))
        PolygonOps.inpolygon(point2d, pts2d) != 0
    end
end

function _above_or_below_polygon(pt::AbstractCoordinatePoint, p::Polygon{N,T}) where {N,T}

    rot = _get_rot_for_rotation_on_xy_plane(p)
    # handle polygons with zero area
    all(isfinite.(rot)) || return false

    vs = vertices(p)

    pts2d = SVector{N+1, SVector{2,T}}(
        SVector{2,T}(view(rot * (vs[i % N + 1] - cartesian_zero), 1:2)) for i in 0:N
    )
    
    rotated_pt = rot * (pt - cartesian_zero)
    point2d = SVector{2,T}(view(rotated_pt, 1:2))
    
    return PolygonOps.inpolygon(point2d, pts2d) != 0
end

function _transform_into_global_coordinate_system(poly::Polygon{N, T}, p::AbstractPrimitive{T}) where {N, T}
    Polygon{N,T}(broadcast(pt -> _transform_into_global_coordinate_system(pt, p), poly.points))
end

function distance(pt::CartesianPoint, p::Polygon)
    return if _above_or_below_polygon(pt, p)
        distance(pt, Plane(p))
    else
        es = edges(p)
        ds = map(e -> distance(pt, e), es)
        minimum(ds)
    end
end

@inline function distance_to_surface(pt::AbstractCoordinatePoint{T}, p::Polygon{N,T})::T where {N,T}
    return distance(CartesianPoint(pt), p)
end
