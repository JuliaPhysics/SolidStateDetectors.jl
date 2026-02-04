"""
    struct ConeMantle{T,TR,TP,D} <: AbstractCurvedSurfacePrimitive{T}

Surface primitive describing the mantle of a [`Cone`](@ref).

## Parametric types
* `T`: Precision type.
* `TR`: Type of the radius `r`.
    * `TR == T`: CylinderMantle (constant radius `r` at all `z`).
    * `TR == Tuple{T, T}`: VaryingCylinderMantle (bottom radius at `r[1]`, top radius at `r[2]`).
* `TP`: Type of the angular range `φ`.
    * `TP == Nothing`: Full 2π Cone.
    * `TP == T`: Partial ConeMantle ranging from `0` to `φ`.
* `D`: Direction in which the normal vector points (`:inwards` or `:outwards`).
    
## Fields
* `r::TR`: Definition of the radius of the `ConeMantle` (in m).
* `φ::TP`: Range in polar angle `φ` over which the `ConeMantle` extends (in radians).
* `hZ::T`: Half of the height of the `ConeMantle` (in m).
* `origin::CartesianPoint{T}`: Origin of the `Cone` which has this `ConeMantle` as surface.
* `rotation::SMatrix{3,3,T,9}`: Rotation matrix of the `Cone` which has this `ConeMantle` as surface.
"""
struct ConeMantle{T,TR,TP<:Union{Nothing,T},D} <: AbstractCurvedSurfacePrimitive{T}
    r::TR
    φ::TP
    hZ::T

    origin::CartesianPoint{T}
    rotation::SMatrix{3,3,T,9}
end

#Type conversion happens here
function ConeMantle{T}(D, r, φ, hZ, origin, rotation) where {T}
    _r = _csg_convert_args(T, r)
    (_φ, _rotation) = _handle_phi(_csg_convert_args(T, φ), rotation)
    _hZ = _csg_convert_args(T, hZ)
    ConeMantle{T,typeof(_r),typeof(_φ),D}(_r, _φ, _hZ, origin, _rotation)
end

#Type promotion happens here
function ConeMantle(D,r::TR, φ::TP, hZ::TZ, origin::PT, rotation::ROT) where {TR, TP, TZ, PT, ROT}
    eltypes = _csg_get_promoted_eltype.((TR, TP, TZ, PT, ROT))
    T = float(promote_type(eltypes...))
    ConeMantle{T}(D, r, φ, hZ, origin, rotation)
end

function ConeMantle(D::Symbol=:outwards;
    # define default parameters as Int to not influence type promotion
    r = 1, 
    φ = nothing,
    hZ = 1,
    origin = zero(CartesianPoint{Int}), 
    rotation = one(SMatrix{3, 3, Int, 9}))
    ConeMantle(D, r, φ, hZ, origin, rotation)
end

function ConeMantle{T}(D::Symbol=:outwards;
    r = 1.0, 
    φ = nothing,
    hZ= 1.0,
    origin = zero(CartesianPoint{Float64}), 
    rotation = one(SMatrix{3, 3, Float64, 9})
) where {T}
    ConeMantle{T}(D, r, φ, hZ, origin, rotation)
end

flip(c::ConeMantle{T,TR,TP,:inwards}) where {T,TR,TP} = 
    ConeMantle{T,TR,TP,:outwards}(c.r, c.φ, c.hZ, c.origin, c.rotation)
flip(c::ConeMantle{T,TR,TP,:outwards}) where {T,TR,TP} = 
    ConeMantle{T,TR,TP,:inwards}(c.r, c.φ, c.hZ, c.origin, c.rotation) 

radius_at_z(hZ::T, rBot::T, rTop::T, z::T) where {T} = iszero(hZ) ? rBot : rBot + (hZ+z)*(rTop - rBot)/(2hZ) 
radius_at_z(cm::ConeMantle{T,T}, z::T) where {T} = cm.r
radius_at_z(cm::ConeMantle{T,Tuple{T,T}}, z::T) where {T} = radius_at_z(cm.hZ, cm.r[1], cm.r[2], z)

get_φ_limits(cm::ConeMantle{T,<:Any,T}) where {T} = T(0), cm.φ
get_φ_limits(cm::ConeMantle{T,<:Any,Nothing}) where {T} = T(0), T(2π)

# function normal(cm::ConeMantle{T,T,<:Any,:inwards}, pt::CartesianPoint{T})::CartesianVector{T} where {T}
#     pto = _transform_into_object_coordinate_system(pt, cm)
#     cyl = CylindricalPoint(pto)
#     return _transform_into_global_coordinate_system(
#             CartesianVector(CartesianPoint(CylindricalPoint{T}(-cyl.r, cyl.φ, zero(T)))), cm)
# end
# function normal(cm::ConeMantle{T,T,<:Any,:outwards}, pt::CartesianPoint{T})::CartesianVector{T} where {T}
#     pto = _transform_into_object_coordinate_system(pt, cm)
#     cyl = CylindricalPoint(pto)
#     return _transform_into_global_coordinate_system(
#             CartesianVector(CartesianPoint(CylindricalPoint{T}(cyl.r, cyl.φ, zero(T)))), cm)
# end

function normal(cm::ConeMantle{T,T,<:Any,:outwards}, pt::CartesianPoint{T})::CartesianVector{T} where {T}
    p_local = _transform_into_object_coordinate_system(pt, cm)
    x = p_local.x
    y = p_local.y
    normal_local = normalize(CartesianVector(x, y, zero(T)))
    return _transform_into_global_coordinate_system(normal_local, cm)
end

function normal(cm::ConeMantle{T,T,<:Any,:inwards}, pt::CartesianPoint{T})::CartesianVector{T} where {T}
    p_local = _transform_into_object_coordinate_system(pt, cm)
    x = p_local.x
    y = p_local.y
    normal_local = -normalize(CartesianVector(x, y, zero(T)))
    return _transform_into_global_coordinate_system(normal_local, cm)
end

# function normal(cm::ConeMantle{T,Tuple{T,T},<:Any,:inwards}, pt::CartesianPoint{T})::CartesianVector{T} where {T}
#     pto = _transform_into_object_coordinate_system(pt, cm)
#     cyl = CylindricalPoint(pto)
#     Δr = cm.r[2] - cm.r[1]
#     Δz = 2cm.hZ
#     return _transform_into_global_coordinate_system(
#             CartesianVector(CartesianPoint(CylindricalPoint{T}(-one(T), cyl.φ, Δr / Δz))), cm)
# end
function normal(cm::ConeMantle{T,Tuple{T,T},<:Any,:inwards}, pt::CartesianPoint{T})::CartesianVector{T} where {T}
    p_local = _transform_into_object_coordinate_system(pt, cm)
    x = p_local.x
    y = p_local.y
    
    Δr = cm.r[2] - cm.r[1]
    Δz = 2cm.hZ
    norm_vec_local = normalize(CartesianVector(-x, -y, Δr / Δz))
    return _transform_into_global_coordinate_system(norm_vec_local, cm)

end

# function normal(cm::ConeMantle{T,Tuple{T,T},<:Any,:outwards}, pt::CartesianPoint{T})::CartesianVector{T} where {T}
#     pto = _transform_into_object_coordinate_system(pt, cm)
#     cyl = CylindricalPoint(pto)
#     Δr = cm.r[2] - cm.r[1]
#     Δz = 2cm.hZ
#     return _transform_into_global_coordinate_system(
#             CartesianVector(CartesianPoint(CylindricalPoint{T}( one(T), cyl.φ, -Δr / Δz))), cm)
# end

function normal(cm::ConeMantle{T,Tuple{T,T},<:Any,:outwards}, pt::CartesianPoint{T})::CartesianVector{T} where {T}
    p_local = _transform_into_object_coordinate_system(pt, cm)
    x = p_local.x
    y = p_local.y

    Δr = cm.r[2] - cm.r[1]
    Δz = 2cm.hZ
    norm_vec_local = normalize(CartesianVector(x, y, -Δr / Δz))
    return _transform_into_global_coordinate_system(norm_vec_local, cm)

end

function vertices(cm::ConeMantle{T}, n_arc::Int)::Vector{CartesianPoint{T}} where {T}
    φMin, φMax = get_φ_limits(cm)
    n_arc = _get_n_points_in_arc_φ(cm, n_arc)
    φ = range(φMin, stop = φMax, length = n_arc+1)
    rbot = radius_at_z(cm,-cm.hZ)
    rtop = radius_at_z(cm,cm.hZ)
    botcircle = [_transform_into_global_coordinate_system(CartesianPoint{T}(rbot*cos(φ), rbot*sin(φ), -cm.hZ), cm) for φ in φ]
    topcircle = [_transform_into_global_coordinate_system(CartesianPoint{T}(rtop*cos(φ), rtop*sin(φ), cm.hZ), cm) for φ in φ]
    append!(botcircle, topcircle)
end

function sample(cm::ConeMantle{T}, spacing::T)::Vector{CartesianPoint{T}} where {T}
    φMin, φMax = get_φ_limits(cm)
    Δφ = abs(φMax - φMin)
    full2π = isnothing(cm.φ)
    rbot = radius_at_z(cm,-cm.hZ)
    rtop = radius_at_z(cm,cm.hZ)
    l = hypot(2cm.hZ, rtop - rbot)
    z = range(-cm.hZ, stop = cm.hZ, length = max(2,1+Int(ceil(l/spacing))))
    [_transform_into_global_coordinate_system(CartesianPoint{T}((radius_at_z(cm,z) .* reverse(sincos(φ)))..., z), cm) 
        for z in z for r in [radius_at_z(cm,z)] 
        for φ in (r == 0 ? [φMin] : range(φMin, stop = φMax - full2π*spacing/radius_at_z(cm,z), length = max(2,1+Int(ceil(Δφ*radius_at_z(cm,z)/spacing)))))]
end

function connections(cm::ConeMantle, n_arc::Int)::Vector{Vector{Int}} 
    n_arc = _get_n_points_in_arc_φ(cm, n_arc)
    [[i,i+1,i+n_arc+2,i+n_arc+1] for i in 1:n_arc]
end

function connections(cm::ConeMantle, n_arc::Int, n_vert_lines::Int)::Vector{Vector{Int}} 
    n_arc = _get_n_points_in_arc_φ(cm, n_arc)
    verts = [[i, i + n_arc + 1] for i in _get_vert_lines_range(cm,n_arc,n_vert_lines)]
    circ1 =  [[i, i + 1] for i in 1:n_arc]
    circ2 =  [[i, i + 1] for i in n_arc+2:2*n_arc+1]
    append!(verts, circ1, circ2)
end

get_label_name(::ConeMantle) = "Cone Mantle"

const FullConeMantle{T,D} = ConeMantle{T,Tuple{T,T},Nothing,D} # ugly name but works for now, should just be `ConeMantle`...
const PartialConeMantle{T,D} = ConeMantle{T,Tuple{T,T},T,D}


extremum(cm::FullConeMantle{T}) where {T} = hypot(cm.hZ, max(cm.r...))

extremum(cm::PartialConeMantle{T}) where {T} = hypot(cm.hZ, max(cm.r...))


function lines(sp::FullConeMantle{T}; n = 2) where {T} 
    bot_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T),zero(T), -sp.hZ), sp)
    bot_ellipse = Circle{T}(sp.r[1], sp.φ, bot_origin, sp.rotation)
    top_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T),zero(T), +sp.hZ), sp)
    top_ellipse = Circle{T}(sp.r[2], sp.φ, top_origin, sp.rotation)
    φs = range(T(0), step = T(2π) / n, length = n)
    edges = [ Edge{T}(_transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[1], φ, -sp.hZ)), sp),
                      _transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[2], φ, +sp.hZ)), sp)) for φ in φs ]
    bot_ellipse, top_ellipse, edges
end
function lines(sp::PartialConeMantle{T}; n = 2) where {T} 
    bot_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -sp.hZ), sp)
    bot_ellipse = PartialCircle{T}(sp.r[1], sp.φ, bot_origin, sp.rotation)
    top_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +sp.hZ), sp)
    top_ellipse = PartialCircle{T}(sp.r[2], sp.φ, top_origin, sp.rotation)
    φs = range(zero(T), stop = sp.φ, length = n)
    edges = [ Edge{T}(_transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[1], φ, -sp.hZ)), sp),
                      _transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[2], φ, +sp.hZ)), sp)) for φ in φs ]
    bot_ellipse, top_ellipse, edges
end


"""
    intersection(cm::ConeMantle{T,Tuple{T,T}}, l::Line{T}) where {T}

Calculates the intersections of a `Line` with a `ConeMantle`.

## Arguments
* `cm::ConeMantle{T,Tuple{T,T}}`: The `ConeMantle`.
* `l::Line{T}`: The `Line`.

!!! note 
    The function will always return 2 CartesianPoint's.
    If the line just touches the mantle, the two points will be the same. 
    If the line does not touch the mantle at all, the two points will have NaN's as there coordinates.
    If the line crosses the mantle only once, two points will be returned. The two points will be the same point (the intersection).
    If the line lies inside the mantle and is parallel to it. The same point will be returned which is the origin of the line. 
"""
function intersection(cm::ConeMantle{T,Tuple{T,T}}, l::Line{T}) where {T}
    obj_l = _transform_into_object_coordinate_system(l, cm) # direction is not normalized
    
    L1 = obj_l.origin.x
    L2 = obj_l.origin.y
    L3 = obj_l.origin.z
    D1 = obj_l.direction.x
    D2 = obj_l.direction.y
    D3 = obj_l.direction.z
    
    hZ = cm.hZ
    R0 = cm.r[1]

    S = (cm.r[2] - cm.r[1]) / 2hZ # r slope

    f1 = D1^2 + D2^2 - D3^2*S^2
    # f1 is zero if there is only one intersection.

    λ1, λ2 = if f1 == 0 # one intersection
        term1 = -hZ^2*S^2 - 2hZ*L3*S^2 - 2*hZ*R0*S + L1^2 + L2^2 - L3^2*S^2 - 2L3*R0*S - R0^2
        term2 = L1*sqrt(D3^2*S^2 - D2^2) - D2*L2 + D3*hZ*S^2 + D3*L3*S^2 + D3*R0*S
        λ = term1 / (2term2) 
        λ, λ
    else # two or no intersections 
        λ = inv(f1)
        term1 = (2D1*L1 + 2D2*L2 - 2D3*hZ*S^2 - 2D3*L3*S^2 - 2D3*R0*S)^2
        term2 = -hZ^2*S^2 - 2hZ*L3*S^2 - 2hZ*R0*S + L1^2 + L2^2 - L3^2*S^2 - 2L3*R0*S - R0^2
        term3 = -D1*L1 - D2*L2 + D3*hZ*S^2 + D3*L3*S^2 + D3*R0*S
        term4 = term1 - 4*f1*term2
        # sq::T = term4 < 0 ? T(NaN) : sqrt(term4) 
        sq::T = sqrt(abs(term4)) 
    
        # @show term4
        λ1 = λ * (-sq/2 + term3) 
        λ2 = λ * (+sq/2 + term3)
        λ1, λ2    
    end
    ints1 = obj_l.origin + λ1 * obj_l.direction 
    ints2 = obj_l.origin + λ2 * obj_l.direction 
    # Could be reintroduced as sanity check
    # ints1 = ints1 * ifelse(abs(ints1.z) <= cm.hZ, one(T), T(NaN))
    # ints2 = ints2 * ifelse(abs(ints2.z) <= cm.hZ, one(T), T(NaN))
    return _transform_into_global_coordinate_system(ints1, cm), 
           _transform_into_global_coordinate_system(ints2, cm)
end


# function get_2d_grid_ticks_and_proj(cm::ConeMantle{T}, t) where {T}
#     pts = extreme_points(cm)
    
#     t_idx_range_x, t_idx_range_y, t_idx_range_z = get_min_max_index_ranges((pts, t))
#     ls = (length(t_idx_range_x), length(t_idx_range_y), length(t_idx_range_z))
#     ls = (
#         ls[1] == 1 ? typemax(eltype(ls)) : ls[1],
#         ls[2] == 1 ? typemax(eltype(ls)) : ls[2],
#         ls[3] == 1 ? typemax(eltype(ls)) : ls[3]
#     )
#     eX = CartesianVector{T}(one(T),zero(T),zero(T))
#     eY = CartesianVector{T}(zero(T),one(T),zero(T))
#     eZ = CartesianVector{T}(zero(T),zero(T),one(T))
#     axis = cm.rotation * CartesianVector{T}(zero(T), zero(T), one(T))

#     proj, t1, t2 = if (axis × eZ != zero(CartesianVector{T}))
#         Val{:xy}(), t_idx_range_x, t_idx_range_y
#     elseif (axis × eY != zero(CartesianVector{T}))
#         Val{:xz}(), t_idx_range_x, t_idx_range_z
#     elseif (axis × eX != zero(CartesianVector{T}))
#         Val{:yz}(), t_idx_range_y, t_idx_range_z
#     else
#         error("Sampling Error. Have to extend cases")
#         # Should never happen. This else case can be removed after testing.
#     end

#     return t1, t2, proj
# end


function distance_to_surface(pt::AbstractCoordinatePoint{T}, c::ConeMantle{T, <:Any, TP, <:Any})::T where {T, TP<:Union{Nothing,T}}

    pt_cart = _transform_into_object_coordinate_system(CartesianPoint(pt), c)
    pt_cyl  = CylindricalPoint(pt_cart)

    rbot::T, rtop::T = c.r isa Tuple ? c.r : (c.r, c.r)
    zMin::T, zMax::T = _linear_endpoints(c.hZ)

    # Generator in the meridional (r–z) plane
    edge_rz = Edge{T}(CartesianPoint{T}(rbot, zero(T), zMin), CartesianPoint{T}(rtop, zero(T), zMax))
    
    # Full cone: rotational symmetry
    return if isnothing(c.φ) || _in_φ(pt_cyl, c.φ)
        distance_to_line(CartesianPoint{T}(pt_cyl.r, zero(T), pt_cyl.z), edge_rz)
    else
        # Nearest boundary: either 0 or c.φ
        φNear = _φNear(pt_cyl.φ, c.φ)
        sinφ, cosφ = sincos(φNear)
        p1 = CartesianPoint{T}(rbot*cosφ, rbot*sinφ, zMin)
        p2 = CartesianPoint{T}(rtop*cosφ, rtop*sinφ, zMax)
        distance_to_line(pt_cart, Edge(p1, p2))
    end
end
