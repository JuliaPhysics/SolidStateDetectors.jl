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
            barycenter(p.points), Plane(p).normal / 5
        end
    end
end