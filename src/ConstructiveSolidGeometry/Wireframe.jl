function get_plot_points(v::AbstractVolumePrimitive{T}) where {T <: AbstractFloat}
    plot_points = Vector{CartesianPoint{T}}[]
    for surf in get_decomposed_surfaces(v)
        append!(plot_points, get_plot_points(surf))
    end
    unique(plot_points)
end
