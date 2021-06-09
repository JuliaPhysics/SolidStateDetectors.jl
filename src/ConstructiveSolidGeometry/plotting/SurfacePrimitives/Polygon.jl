@recipe function f(p::Polygon; show_normal = false, normallength = 0, normalcolor = missing)
    @series begin
        x = map(p -> p[1], [p.points..., p.points[1]])
        y = map(p -> p[2], [p.points..., p.points[1]])
        z = map(p -> p[3], [p.points..., p.points[1]])
        x, y, z
    end
    if show_normal
        @series begin
            seriestype := :vector
            if !ismissing(normalcolor) 
                linecolor := normalcolor
            end
            mean(p.points), Plane(p).normal * normallength
        end
    end
end

@recipe function f(vp::AbstractVector{<:Polygon})
    for p in vp
        @series begin
            show_normal := true
            p
        end
    end
end

