@recipe function f(c::AbstractGeometry{T}) where {T}
    
    seriescolor --> :orange
    linewidth --> 2
    
    @series begin
        label --> "Geometry"
        []
    end

    for ls in get_plot_points(c)
        @series begin
            label := ""
            ls
        end
    end
end

@recipe function f(ls::LineSegments{T}) where {T}
    points::Vector{CartesianPoint{T}} = ls.points
    map(p -> p.x, points), map(p -> p.y, points), map(p -> p.z, points)
end


