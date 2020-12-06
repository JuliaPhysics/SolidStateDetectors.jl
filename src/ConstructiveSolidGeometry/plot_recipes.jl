@recipe function f(g::AbstractGeometry{T}) where {T}
    
    seriescolor --> :orange
    linewidth --> 2
    
    @series begin
        label --> "Geometry"
        []
    end
    
    pos, _ = get_decomposed_volumes(g)
    
    for p in pos
        for points in get_plot_points(p)
            @series begin
                label := ""
                points
            end
        end
    end
end


@recipe function f(points::Vector{CartesianPoint{T}}) where {T}
    map(p -> p.x, points), map(p -> p.y, points), map(p -> p.z, points)
end


