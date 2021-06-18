abstract type AbstractPlanarSurfacePrimitive{T} <: AbstractSurfacePrimitive{T} end
abstract type AbstractCurvedSurfacePrimitive{T} <: AbstractSurfacePrimitive{T} end

include("Plane.jl")
include("Polygon.jl")
include("ConeMantle.jl")
include("EllipticalSurface.jl")


function get_min_max_index_ranges(a::Union{
        <:Tuple{AbstractVector{<:CartesianPoint{T}}, CartesianTicksTuple{T}},
        <:Tuple{AbstractVector{<:CylindricalPoint{T}}, CylindricalTicksTuple{T}}
        }) where {T}
    pts, t = a[1], a[2]
    ax1, ax2, ax3 = getindex.(pts, 1), getindex.(pts, 2), getindex.(pts, 3)
    ax1_min_idx = searchsortedfirst(t[1], minimum(ax1))
    ax2_min_idx = searchsortedfirst(t[2], minimum(ax2))
    ax3_min_idx = searchsortedfirst(t[3], minimum(ax3))
    ax1_max_idx = searchsortedfirst(t[1], maximum(ax1))
    ax2_max_idx = searchsortedfirst(t[2], maximum(ax2))
    ax3_max_idx = searchsortedfirst(t[3], maximum(ax3))
    ls = (length(t[1]), length(t[2]), length(t[3]))
    if ax1_max_idx > ls[1] ax1_max_idx = ls[1] end
    if ax2_max_idx > ls[2] ax2_max_idx = ls[2] end
    if ax3_max_idx > ls[3] ax3_max_idx = ls[3] end
    t_idx_range_ax1 = ax1_min_idx:ax1_max_idx
    t_idx_range_ax2 = ax2_min_idx:ax2_max_idx
    t_idx_range_ax3 = ax3_min_idx:ax3_max_idx
    t_idx_range_ax1, t_idx_range_ax2, t_idx_range_ax3
end

"""
    get_2d_grid_ticks_and_proj(p::AbstractPlanarSurfacePrimitive, t::CartesianTicksTuple{T}) where {N, T}

This function determines the two best dimensions to sample/paint the surface p. 
E.g. `x` & `y` -> `proj = Val{:xy}()`.
The dimensions are picked such that the number of points to evaluate is minimal. 
However, the polygon is not allowed to be parallel to the remaining dimension, e.g. `z`,
because, than, there would not be a single value, but infinite ones,
for `z` in the evalution.
"""
function get_2d_grid_ticks_and_proj(p::AbstractPlanarSurfacePrimitive, t::CartesianTicksTuple{T}) where {N, T}
    # This method would actually work for any flat surface, e.g. elipse
    pts = extreme_points(p)
    t_idx_range_x, t_idx_range_y, t_idx_range_z = get_min_max_index_ranges((pts, t))
    ls = (length(t_idx_range_x), length(t_idx_range_y), length(t_idx_range_z))
    ls = (
        ls[1] == 1 ? typemax(eltype(ls)) : ls[1],
        ls[2] == 1 ? typemax(eltype(ls)) : ls[2],
        ls[3] == 1 ? typemax(eltype(ls)) : ls[3]
    )
    n = normal(p)
    proj, t1, t2 = if (n ⋅ CartesianVector{T}(zero(T),zero(T),one(T)) != 0)
        Val{:xy}(), t_idx_range_x, t_idx_range_y
    elseif (n ⋅ CartesianVector{T}(zero(T),one(T),zero(T)) != 0)
        Val{:xz}(), t_idx_range_x, t_idx_range_z
    elseif n ⋅ CartesianVector{T}(one(T),zero(T), zero(T)) != 0
        Val{:yz}(), t_idx_range_y, t_idx_range_z
    else
        error("Sampling Error. Have to extend cases")
        # Should never happen. This else case can be removed after testing.
    end
    return t1, t2, proj
end


"""
    get_2d_grid_ticks_and_proj(p::AbstractPlanarSurfacePrimitive, t::CylindricalTicksTuple{T}) where {N, T}

This function determines the two best dimensions to sample/paint the surface p. 
E.g. `r` & `φ` -> `proj = Val{:rφ}()`.
The dimensions are picked such that the number of points to evaluate is minimal. 
However, the polygon is not allowed to be parallel to the remaining dimension, e.g. `z`,
because, than, there would not be a single value, but infinite ones, for `z` in the evalution.
"""
function get_2d_grid_ticks_and_proj(p::AbstractPlanarSurfacePrimitive, t::CylindricalTicksTuple{T}) where {N, T}
    # This method would actually work for any flat surface, e.g. elipse
    pts = extreme_points(p)
    t_idx_range_r, t_idx_range_φ, t_idx_range_z = get_min_max_index_ranges((CylindricalPoint.(pts), t))
    ls = (length(t_idx_range_r), length(t_idx_range_φ), length(t_idx_range_z))
    ls = (
        ls[1] == 1 ? typemax(eltype(ls)) : ls[1],
        ls[2] == 1 ? typemax(eltype(ls)) : ls[2],
        ls[3] == 1 ? typemax(eltype(ls)) : ls[3]
    )
    n = normal(p)
    p_r = pts[findmax(broadcast(p -> abs(p[1]), pts))[2]]
    # We skip the `:rz`-case as there could be two intersections with the arc and the polygon.
    proj, t1, t2 = if (n ⋅ CartesianVector{T}(zero(T),zero(T),one(T)) != 0)
        # evaluate the z value -> same as for the cartesian case
        Val{:rφ}(), t_idx_range_r, t_idx_range_φ
    elseif (n ⋅ p_r != 0)
        Val{:φz}(), t_idx_range_φ, t_idx_range_z
    else
        Val{:rz}(), t_idx_range_r, t_idx_range_z
    end
    return t1, t2, proj
end
