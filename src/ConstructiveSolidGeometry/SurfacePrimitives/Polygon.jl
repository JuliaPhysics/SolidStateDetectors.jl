
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

function get_2d_grid_ticks_and_proj(p::Polygon{N, T}, t::CartesianTicksTuple{T}) where {N, T}
    # This method would actually work for any flat surface, e.g. elipse
    xs, ys, zs = getindex.(p.points, 1), getindex.(p.points, 2), getindex.(p.points, 3)
    xmin_idx = searchsortedfirst(t[1], minimum(xs))
    ymin_idx = searchsortedfirst(t[2], minimum(ys))
    zmin_idx = searchsortedfirst(t[3], minimum(zs))
    xmax_idx = searchsortedfirst(t[1], maximum(xs))
    ymax_idx = searchsortedfirst(t[2], maximum(ys))
    zmax_idx = searchsortedfirst(t[3], maximum(zs))
    ls = (length(t[1]), length(t[2]), length(t[3]))
    if xmax_idx > ls[1] zmax_idx = ls[1] end
    if ymax_idx > ls[2] zmax_idx = ls[2] end
    if zmax_idx > ls[3] zmax_idx = ls[3] end
    t_idx_range_x = xmin_idx:xmax_idx
    t_idx_range_y = ymin_idx:ymax_idx
    t_idx_range_z = zmin_idx:zmax_idx
    ls = (length(t_idx_range_x), length(t_idx_range_y), length(t_idx_range_z))
    ls = (
        ls[1] == 1 ? typemax(eltype(ls)) : ls[1],
        ls[2] == 1 ? typemax(eltype(ls)) : ls[2],
        ls[3] == 1 ? typemax(eltype(ls)) : ls[3]
    )
    n = normal(p)
    proj, t1, t2 = if ls[1] < ls[3] && ls[2] < ls[3] && 
                      (n ⋅ CartesianVector{T}(zero(T),zero(T),one(T)) != 0)
        Val{:xy}(), t_idx_range_x, t_idx_range_y
    elseif ls[1] < ls[2] && ls[3] < ls[2] && 
                    (n ⋅ CartesianVector{T}(zero(T),one(T),zero(T)) != 0)
        Val{:xz}(), t_idx_range_x, t_idx_range_z
    elseif n ⋅ CartesianVector{T}(one(T),zero(T), zero(T)) != 0
        Val{:yz}(), t_idx_range_y, t_idx_range_z
    else
        error("Sampling Error. Have to extend cases")
    end
    return t1, t2, proj
end
function get_2d_grid_ticks_and_proj(p::Polygon{N, T}, t::CylindricalTicksTuple{T}) where {N, T}
    # This method would actually work for any flat surface, e.g. elipse
    cyls = CylindricalPoint.(p.points)
    rs, φs, zs = getindex.(cyls, 1), getindex.(cyls, 2), getindex.(cyls, 3)
    rmin_idx = searchsortedfirst(t[1], minimum(rs))
    φmin_idx = searchsortedfirst(t[2], minimum(φs))
    zmin_idx = searchsortedfirst(t[3], minimum(zs))
    rmax_idx = searchsortedfirst(t[1], maximum(rs))
    φmax_idx = searchsortedfirst(t[2], maximum(φs))
    zmax_idx = searchsortedfirst(t[3], maximum(zs))
    ls = (length(t[1]), length(t[2]), length(t[3]))
    if rmax_idx > ls[1] rmax_idx = ls[1] end
    if φmax_idx > ls[2] φmax_idx = ls[2] end
    if zmax_idx > ls[3] zmax_idx = ls[3] end
    t_idx_range_r = rmin_idx:rmax_idx
    t_idx_range_φ = φmin_idx:φmax_idx
    t_idx_range_z = zmin_idx:zmax_idx
    ls = (length(t_idx_range_r), length(t_idx_range_φ), length(t_idx_range_z))
    ls = (
        ls[1] == 1 ? typemax(eltype(ls)) : ls[1],
        ls[2] == 1 ? typemax(eltype(ls)) : ls[2],
        ls[3] == 1 ? typemax(eltype(ls)) : ls[3]
    )
    n = normal(p)
    proj, t1, t2 = if ls[1] < ls[3] && ls[2] < ls[3] && 
                      (n ⋅ CartesianVector{T}(zero(T),zero(T),one(T)) != 0)
        Val{:rφ}(), t_idx_range_r, t_idx_range_φ
    # elseif ls[1] < ls[2] && ls[3] < ls[2] # for a polygon we do not need Val{:rz}
    #     Val{:rz}(), t_idx_range_r, t_idx_range_z
    else
        Val{:φz}(), t_idx_range_φ, t_idx_range_z
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

