function _get_tessellation_nodes(nodes::Vector{SVector{2,T}}) where {T}
    u = map(p -> p[1], nodes)
    min_u = minimum(u)
    width_u = maximum(u) - min_u
    v = map(p -> p[2], nodes)
    min_v = minimum(v)
    width_v = maximum(v) - min_v
    width = VoronoiDelaunay.max_coord - VoronoiDelaunay.min_coord
    #max_coord and min_coord are defined and exported in VoronoiDelaunay package: 1 + eps(Float64), 2 - 2eps(Float64).
    Vector{Point2D}(map(p -> VoronoiDelaunay.Point2D(width*(p[1] - min_u)/width_u + min_coord, width*(p[2] - min_v)/width_v + min_coord), nodes)), min_u, width_u, min_v, width_v, width
end

_get_planar_point_from_tesselation_node(node::VoronoiDelaunay.Point2D, min_u::Real, width_u::Real, min_v::Real, width_v::Real, width::Real, T::Type) = SVector{2,T}(width_u*(node._x - min_coord)/width + min_u, width_v*(node._y - min_coord)/width + min_v)

function triangles(p::Polygon{N,T})::Vector{Polygon{3,T}} where {N, T}
    rot = _get_rot_for_rotation_on_xy_plane(p)
    invrot = inv(rot)
    vs = vertices(p)
    pts2d = SVector{N+1, SVector{2, T}}(
        view(rot * vs[i%N + 1], 1:2) for i in 0:N  
    )
    tess = DelaunayTessellation(N)
    tess_nodes, min_u, width_u, min_v, width_v, width = _get_tessellation_nodes(pts2d[1:end-1])
    push!(tess, tess_nodes)
    triangles = Vector{Polygon{3,T}}()
    for dt in tess
        tri = SVector{3, SVector{2, T}}(_get_planar_point_from_tesselation_node(v, min_u, width_u, min_v, width_v, width, T) for v in (dt._a, dt._b, dt._c))
        if PolygonOps.inpolygon(sum(tri)/3, pts2d) != 0
            tri3d = SVector{3, CartesianPoint{T}}(CartesianPoint{T}(invrot*[p[1],p[2],0]) for p in tri)
            if !isapprox(cross(tri[3]-tri[1], tri[2]-tri[1]), 0, atol = 0.00000001)
                push!(triangles, Polygon{3,T}(tri3d))
            end
        end
    end
    triangles
end

function trimesh(tri::Polygon{3,T})::PolyMesh{T,3} where {T}
    if isapprox(normal(tri) ⋅ [0,0,1], 0, atol = 0.001)               
        tri = _rotate_triangle(tri,RotXY(0.0000001,0.0000001))
    end
    p1, p2, p3 = vertices(tri)  
    PolyMesh{T,3}([[p1.x, p2.x, p3.x]], [[p1.y, p2.y, p3.y]], [[p1.z, p2.z, p3.z]])
end

function _rotate_triangle(tri::Polygon{3,T}, R::AbstractMatrix)::Polygon{3,T} where {T}
    vs = vertices(tri)
    Polygon{3,T}([R*vs[1],R*vs[2],R*vs[3]])
end

function trimesh(p::Polygon{N,T})::PolyMesh{T,3} where {N, T}
    tris = triangles(p)
    if isapprox(normal(p) ⋅ [0,0,1], 0, atol = 0.001)               
        tris = [_rotate_triangle(t,RotXY(0.0000001,0.0000001)) for t in tris]
    end    
    x = [SVector{3,T}(t.points[1].x, t.points[2].x, t.points[3].x) for t in tris]
    y = [SVector{3,T}(t.points[1].y, t.points[2].y, t.points[3].y) for t in tris]
    z = [SVector{3,T}(t.points[1].z, t.points[2].z, t.points[3].z) for t in tris]
    PolyMesh{T,3}(x, y, z)
end

@recipe function f(p::Polygon)
    colorbar := false
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :surface
        @series begin
            label --> "Polygon"
            trimesh(p)
        end
    else
        linecolor --> :black
        @series begin
            label --> "Polygon"
            x = map(p -> p[1], [p.points..., p.points[1]])
            y = map(p -> p[2], [p.points..., p.points[1]])
            z = map(p -> p[3], [p.points..., p.points[1]])
            x, y, z
        end
    end
    if !haskey(plotattributes, :show_normal) || plotattributes[:show_normal]
        @series begin
            label := nothing
            seriestype := :vector
            mean(p.points), Plane(p).normal / 5
        end
    end
end