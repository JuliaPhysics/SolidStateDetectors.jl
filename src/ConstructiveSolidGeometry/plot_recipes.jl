@recipe function f(c::AbstractGeometry{T}) where {T}
    
    seriescolor --> :orange
    linewidth --> 2
    
    @series begin
        label --> "Geometry"
        []
    end

    for points in get_plot_points(c)
        @series begin
            label := ""
            points
        end
    end
end

@recipe function f(points::Vector{CartesianPoint{T}}) where {T}
    map(p -> p.x, points), map(p -> p.y, points), map(p -> p.z, points)
end


