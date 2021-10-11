function mesh(p::Polygon{N,T})::Mesh{T} where {N,T}
    #assumes convex polygonsgr()
    vs = vertices(p)
    x, y, z = broadcast(i -> getindex.(vs, i), (1,2,3))
    connections = [collect(1:N)]
    Mesh{T}(x,y,z,connections)
end

@recipe function f(p::Polygon)
    seriestype --> :mesh3d
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :mesh3d
        @series begin
            label --> "Polygon"
            mesh(p)
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