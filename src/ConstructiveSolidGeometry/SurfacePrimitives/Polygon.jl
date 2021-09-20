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

extreme_points(p::Polygon) = p.points

function edges(p::Triangle{T}) where {T}
    vs = vertices(p)
    return (
        Edge(vs[1], vs[2]),
        Edge(vs[2], vs[3]),
        Edge(vs[3], vs[1])
    )
end
function edges(p::Quadrangle{T}) where {T}
    vs = vertices(p)
    return (
        Edge(vs[1], vs[2]),
        Edge(vs[2], vs[3]),
        Edge(vs[3], vs[4]),
        Edge(vs[4], vs[1])
    )
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

function in(pt::CartesianPoint{T}, p::Polygon{N,T}) where {N,T}
    plane = Plane(p)
    if (pt - origin(plane)) ⋅ normal(plane) == 0
        rot = _get_rot_for_rotation_on_xy_plane(p)
        vs = vertices(p)
        pts2d = SVector{N+1, SVector{2, T}}(
            view(rot * vs[i%N + 1], 1:2) for i in 0:N  
        )
        # PolygonOps.inpolygon -> in = 1, on = -1, out = 0)
        b = PolygonOps.inpolygon(view(rot * pt, 1:2), pts2d) != 0 
    end
    return b
end

function _above_or_below_polygon(pt::AbstractCoordinatePoint, p::Quadrangle{T}) where {T}
    rot = _get_rot_for_rotation_on_xy_plane(p)
    vs = vertices(p)
    pts2d = SVector{5, SVector{2, T}}(
        view(rot * vs[1], 1:2), 
        view(rot * vs[2], 1:2), 
        view(rot * vs[3], 1:2), 
        view(rot * vs[4], 1:2), 
        view(rot * vs[1], 1:2), 
    )
    # PolygonOps.inpolygon -> in = 1, on = -1, out = 0)
    return PolygonOps.inpolygon(view(rot * pt, 1:2), pts2d) != 0 
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
