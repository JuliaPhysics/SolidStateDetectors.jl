include("PointsAndVectors.jl")

include("VolumePrimitives/VolumePrimitives.jl")

@recipe function f(g::AbstractGeometry{T}; SSD_style = :wireframe, n = 30) where {T}

    seriescolor --> :orange

    @series begin
        label --> "Geometry"
        []
    end

    pos, _ = get_decomposed_volumes(g)

    for p in pos
        if SSD_style == :wireframe
            linewidth --> 2
            for points in get_plot_points(p, n = n)
                @series begin
                    label := ""
                    points
                end
            end
        elseif SSD_style == :surface
            linewidth --> 0.1
            for mesh in get_plot_meshes(p, n = n)
                @series begin
                    label := ""
                    mesh
                end
            end
        end
    end
end


@recipe function f(points::Vector{CartesianPoint{T}}) where {T}
    map(p -> p.x, points), map(p -> p.y, points), map(p -> p.z, points)
end

@recipe function f(points::Vector{CylindricalPoint{T}}) where {T}
    points = CartesianPoint.(points)
    map(p -> p.x, points), map(p -> p.y, points), map(p -> p.z, points)
end

@recipe function f(m::Mesh{T}) where {T}
    seriestype := :surface
    linewidth --> 0.1
    linecolor --> :white
    seriescolor --> :blue
    colorbar := false
    m.x, m.y, m.z
end
