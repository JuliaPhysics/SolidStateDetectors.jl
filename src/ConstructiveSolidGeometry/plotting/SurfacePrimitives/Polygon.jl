@recipe function f(p::Polygon; normallength = 0)
    linecolor --> :black
    @series begin
        label --> "Polygon"
        x = map(p -> p[1], [p.points..., p.points[1]])
        y = map(p -> p[2], [p.points..., p.points[1]])
        z = map(p -> p[3], [p.points..., p.points[1]])
        x, y, z
    end
    if normallength > 0
        @series begin
            label --> "Normal"
            seriestype := :vector
            mean(p.points), Plane(p).normal * normallength
        end
    end
end



