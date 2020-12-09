@recipe function f(c::AbstractGeometry{T}; SSD_style = :wireframe) where {T}

    seriescolor --> :orange

    @series begin
        label --> "Geometry"
        []
    end

    if SSD_style == :wireframe
        linewidth --> 2
        for points in get_plot_points(c)
            @series begin
                label := ""
                points
            end
        end
    elseif SSD_style == :surface
        linewidth --> 0.1
        for mesh in get_plot_meshes(c)
            @series begin
                label := ""
                mesh
            end
        end
    end
end

@recipe function f(points::Vector{CartesianPoint{T}}) where {T}
    map(p -> p.x, points), map(p -> p.y, points), map(p -> p.z, points)
end
