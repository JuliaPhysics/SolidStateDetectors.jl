function _rotate_triangle(tri::Triangle{T}, R::AbstractMatrix)::Triangle{T} where {T}
    vs = vertices(tri)
    Polygon{3,T}([R*vs[1],R*vs[2],R*vs[3]])
end

function triangles(p::Polygon{N,T})::SVector{N-2, Triangle{T}} where {N, T}
    #assumes convex polygons
    vs = vertices(p)
    SVector{N-2, Triangle{T}}(Triangle{T}([vs[1], vs[i+1], vs[i+2]]) for i in 1:N-2)
end

triangles(tri::Triangle{T}) where {T} = SVector{1, Triangle{T}}(tri)

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
    if haskey(plotattributes, :show_normal) && plotattributes[:show_normal]
        @series begin
            label := nothing
            seriestype := :vector
            mean(p.points), Plane(p).normal / 5
        end
    end
end