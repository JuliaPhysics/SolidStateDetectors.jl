function get_sub_ind_ranges(a::Union{
        <:Tuple{AbstractVector{<:CartesianPoint{T}}, CartesianTicksTuple{T}},
        <:Tuple{AbstractVector{<:CylindricalPoint{T}}, CylindricalTicksTuple{T}}
        }) where {T}
    # For Cylindrical Coordinates this only works for z right now
    # There, we will only use t_idx_range_ax3 in the painting (for now...)
    pts, t = a[1], a[2]
    ax1, ax2, ax3 = getindex.(pts, 1), getindex.(pts, 2), getindex.(pts, 3)
    # ± 2 because due to merging close ticks of an axis of a grid,
    # the grid points which should "belong" to corners of a primitive
    # might be a bit outside of the intervals defined by the points 
    # of the `extreme_points`-method.
    # ± 1 might also be enough -> testing...
    ax1_min_idx = searchsortedfirst(t[1], minimum(ax1)) - 2
    ax2_min_idx = searchsortedfirst(t[2], minimum(ax2)) - 2
    ax3_min_idx = searchsortedfirst(t[3], minimum(ax3)) - 2
    ax1_max_idx = searchsortedfirst(t[1], maximum(ax1)) + 2 
    ax2_max_idx = searchsortedfirst(t[2], maximum(ax2)) + 2 
    ax3_max_idx = searchsortedfirst(t[3], maximum(ax3)) + 2 
    ls = (length(t[1]), length(t[2]), length(t[3]))
    if ax1_min_idx < 1 ax1_min_idx = 1 end
    if ax2_min_idx < 1 ax2_min_idx = 1 end
    if ax3_min_idx < 1 ax3_min_idx = 1 end
    if ax1_max_idx > ls[1] ax1_max_idx = ls[1] end
    if ax2_max_idx > ls[2] ax2_max_idx = ls[2] end
    if ax3_max_idx > ls[3] ax3_max_idx = ls[3] end
    t_idx_range_ax1 = ax1_min_idx:ax1_max_idx
    t_idx_range_ax2 = ax2_min_idx:ax2_max_idx
    t_idx_range_ax3 = ax3_min_idx:ax3_max_idx
    t_idx_range_ax1, t_idx_range_ax2, t_idx_range_ax3
end

get_sub_ind_ranges(p::ConstructiveSolidGeometry.AbstractSurfacePrimitive{T}, grid::CartesianGrid{T}) where {N,T} = 
    get_sub_ind_ranges((ConstructiveSolidGeometry.extreme_points(p), TicksTuple(grid)))
get_sub_ind_ranges(p::ConstructiveSolidGeometry.AbstractSurfacePrimitive{T}, grid::CylindricalGrid{T}) where {N,T} = 
    get_sub_ind_ranges((CylindricalPoint.(ConstructiveSolidGeometry.extreme_points(p)), TicksTuple(grid)))

include("Curved.jl")