
struct Polygon{N,T} <: AbstractPlane{T}
    points::SVector{N, CartesianPoint{T}} 
end
const Triangle{T} = Polygon{3, T}
const Quadrangle{T} = Polygon{4, T}

normal(p::Polygon) = normalize((p.points[2] - p.points[1]) × (p.points[3] - p.points[1]))

vertices(p::Polygon) = p.points


function edges(p::Triangle{T}) where {T}
    vs = vertices(p)
    return SVector{3,Edge{T}}(
        Edge(vs[1], vs[2]),
        Edge(vs[2], vs[3]),
        Edge(vs[3], vs[1])
    )
end
function edges(p::Quadrangle{T}) where {T}
    vs = vertices(p)
    return SVector{4,Edge{T}}(
        Edge(vs[1], vs[2]),
        Edge(vs[2], vs[3]),
        Edge(vs[3], vs[4]),
        Edge(vs[4], vs[1])
    )
end

Plane(p::Polygon{N, T}) where {N, T} = Plane{T}(p.points[1], (p.points[2] - p.points[1]) × (p.points[3] - p.points[1]))

function _get_rot_for_rotation_on_xy_plane(p::Polygon{<:Any, T}) where {T}
    n = normal(p)
    rot_angle_rad = CartesianVector{T}(0,0,1) ⋅ n
    return if abs(rot_angle_rad) == 1
        RotationVec(zero(T), zero(T), zero(T))
    else
        rot_angle = -acos(rot_angle_rad)
        rot_axis = rot_angle * normalize(CartesianVector{T}(0,0,1) × n)
        RotationVec(rot_axis[1], rot_axis[2], rot_axis[3])
    end
end

function _rotate_on_xy_plane(p::Polygon{<:Any, T}) where {T}
    rot = _get_rot_for_rotation_on_xy_plane(p)
    map(pt -> rot * pt, p.points)
end


function in(pt::CartesianPoint{T}, p::Quadrangle{T}) where {T}
    b::Bool = in(pt, Plane(p)) 
    if b
        rot = _get_rot_for_rotation_on_xy_plane(p)
        vs = vertices(p)
        pts2d = SVector{5, SVector{2, T}}(
            view(rot * vs[1], 1:2), 
            view(rot * vs[2], 1:2), 
            view(rot * vs[3], 1:2), 
            view(rot * vs[4], 1:2), 
            view(rot * vs[1], 1:2), 
        )
        # PolygonOps.inpolygon -> in = 1, on = -1, out = 0)
        b = PolygonOps.inpolygon(view(rot * pt, 1:2), pts2d) != 0 
    end
    return b
end

function _above_or_below_polygon(pt::AbstractCoordinatePoint, p::Quadrangle{T}) where {T}
    rot = _get_rot_for_rotation_on_xy_plane(p)
    vs = vertices(p)
    pts2d = SVector{5, SVector{2, T}}(
        view(rot * vs[1], 1:2), 
        view(rot * vs[2], 1:2), 
        view(rot * vs[3], 1:2), 
        view(rot * vs[4], 1:2), 
        view(rot * vs[1], 1:2), 
    )
    # PolygonOps.inpolygon -> in = 1, on = -1, out = 0)
    return PolygonOps.inpolygon(view(rot * pt, 1:2), pts2d) != 0 
end

function distance(pt::CartesianPoint, p::Polygon)
    return if _above_or_below_polygon(pt, p)
        distance(pt, Plane(p))
    else
        es = edges(p)
        ds = map(e -> distance(pt, e), es)
        minimum(ds)
    end
end

function _filter_points(pts::AbstractArray, p::Quadrangle{T}) where {T}
    # this functions assume that the points, pts, all lie in the plane of the Quadrangle, p
    rot = _get_rot_for_rotation_on_xy_plane(p)
    vs = vertices(p)
    pts2d = SVector{5, SVector{2, T}}(
        view(rot * vs[1], 1:2), 
        view(rot * vs[2], 1:2), 
        view(rot * vs[3], 1:2), 
        view(rot * vs[4], 1:2), 
        view(rot * vs[1], 1:2), 
    )
    # PolygonOps.inpolygon -> in = 1, on = -1, out = 0)
    filter(pt -> PolygonOps.inpolygon(view(rot * pt, 1:2), pts2d) != 0, pts)
end

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
    get_2d_grid_ticks_and_proj(p::Polygon{N, T}, t::CartesianTicksTuple{T}) where {N, T}

This function determines the two best dimensions to sample/paint the surface. 
E.g. `x` & `y` -> `proj = Val{:xy}()`.
The dimensions are picked such that the number of points to evaluate is minimal. 
However, the polygon is not allowed to be parallel to the remaining dimension, e.g. `z`,
because, than, there would not be a single value, but infinite ones,
for `z` in the evalution.
"""
function get_2d_grid_ticks_and_proj(p::Polygon{N, T}, t::CartesianTicksTuple{T}) where {N, T}
    # This method would actually work for any flat surface, e.g. elipse
    t_idx_range_x, t_idx_range_y, t_idx_range_z = get_min_max_index_ranges((p.points, t))
    ls = (length(t_idx_range_x), length(t_idx_range_y), length(t_idx_range_z))
    ls = (
        ls[1] == 1 ? typemax(eltype(ls)) : ls[1],
        ls[2] == 1 ? typemax(eltype(ls)) : ls[2],
        ls[3] == 1 ? typemax(eltype(ls)) : ls[3]
    )
    n = normal(p)
    proj, t1, t2 = if ls[1] < ls[3] && ls[2] < ls[3] && (n ⋅ CartesianVector{T}(zero(T),zero(T),one(T)) != 0)
        Val{:xy}(), t_idx_range_x, t_idx_range_y
    elseif ls[1] < ls[2] && ls[3] < ls[2] && (n ⋅ CartesianVector{T}(zero(T),one(T),zero(T)) != 0)
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
    get_2d_grid_ticks_and_proj(p::Polygon{N, T}, t::CylindricalTicksTuple{T}) where {N, T}

This function determines the two best dimensions to sample/paint the surface. 
E.g. `r` & `φ` -> `proj = Val{:rφ}()`.
The dimensions are picked such that the number of points to evaluate is minimal. 
However, the polygon is not allowed to be parallel to the remaining dimension, e.g. `z`,
because, than, there would not be a single value, but infinite ones, for `z` in the evalution.
The Val{:rz}-case is excluded as there are two crossections with the arc-line and the polygon. 
"""
function get_2d_grid_ticks_and_proj(p::Polygon{N, T}, t::CylindricalTicksTuple{T}) where {N, T}
    # This method would actually work for any flat surface, e.g. elipse
    t_idx_range_r, t_idx_range_φ, t_idx_range_z = get_min_max_index_ranges((CylindricalPoint.(p.points), t))
    ls = (length(t_idx_range_r), length(t_idx_range_φ), length(t_idx_range_z))
    ls = (
        ls[1] == 1 ? typemax(eltype(ls)) : ls[1],
        ls[2] == 1 ? typemax(eltype(ls)) : ls[2],
        ls[3] == 1 ? typemax(eltype(ls)) : ls[3]
    )
    n = normal(p)
    n_cyl = CylindricalVector{T}(CylindricalPoint(CartesianPoint{T}(n)))
    proj, t1, t2 = if ls[1] < ls[3] && ls[2] < ls[3] && (n ⋅ CartesianVector{T}(zero(T),zero(T),one(T)) != 0)
        # evaluate the z value -> same as for the cartesian case
        Val{:rφ}(), t_idx_range_r, t_idx_range_φ
    elseif n_cyl ⋅ CylindricalVector{T}(one(T),zero(T),zero(T)) != 0
        # evaluate the r value 
        Val{:φz}(), t_idx_range_φ, t_idx_range_z
    else#if ls[1] < ls[2] && ls[3] < ls[2] && (n ⋅ CylindricalVector{T}(zero(T),one(T),zero(T)) != 0)
        # evaluate the φ value -> Issue, there are two solution for φ 
        # Thus, evaluate(p, r, z, Val{:rz}()) will (or might?) return two points. 
        # Therefore we try so skip this by not demanding that 
        # ls[2] < ls[1] && ls[3] < ls[1] for the Val{:φz}()-case.  
        error("Sampling Error. Have to extend cases")
        Val{:rz}(), t_idx_range_r, t_idx_range_z
    end
    return t1, t2, proj
end




function sample(p::Polygon{N, T}, t::CartesianTicksTuple{T}) where {N, T}
    plane = Plane(p)
    t1, t2, proj  = get_2d_grid_ticks_and_proj(p, t)
    samples = Array{CartesianPoint{T}}(undef, length(t1), length(t2))
    for j in eachindex(t2)
        for i in eachindex(t1)
            samples[i, j] = evaluate(plane, t1[i], t2[j], proj)
        end
    end
    _filter_points(Base.ReshapedArray(samples, (length(samples),), ()), p)
end

