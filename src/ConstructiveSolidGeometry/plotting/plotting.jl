include("PointsAndVectors.jl")

include("Meshing.jl")

include("Wireframe.jl")

include("SurfacePrimitives/SurfacePrimitives.jl")

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
            linewidth --> 0.2
            for mesh in get_plot_meshes(p, n = n)
                @series begin
                    label := ""
                    mesh
                end
            end
        end
    end
end

@recipe function f(g::AbstractConstructiveGeometry{T}; SSD_style = :wireframe, n = 30) where {T}
    seriescolor --> :orange

    @series begin
        label --> "Geometry"
        []
    end

    if SSD_style == :wireframe #update to only plot real surfaces
        pos, _ = get_decomposed_volumes(g)
        for p in pos
            linewidth --> 2
            for points in get_plot_points(p, n = n)
                @series begin
                    label := ""
                    points
                end
            end
        end
    elseif SSD_style == :surface
        linewidth --> 0.2
        for mesh in get_plot_meshes(g, n = n)
            @series begin
                label := ""
                mesh
            end
        end
    end
end

@recipe function f(ga::Array{<:AbstractSurfacePrimitive}; n = 30, seriescolor = missing)
    ccolor = ismissing(seriescolor) ? [:blue] : seriescolor
    if !(typeof(ccolor) <: AbstractArray) ccolor = [ccolor] end
    for (cn,g) in enumerate(ga)
        @series begin
            seriescolor := ccolor[(cn-1)%length(ccolor)+1]
            label := ""
            mesh(g, n = n)
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
