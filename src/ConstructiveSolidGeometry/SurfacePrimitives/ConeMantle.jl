"""
    struct ConeMantle{T,TR,TP,D} <: AbstractSurfacePrimitive{T}

T: Type of values, e.g. Float64

* `r::TR`: 
    * TR = Real -> Cylinder Mantle (a = b = r)
    * TR = (Real, Real) -> Cone Mantle (r_bot = r[1], r_top = r[2]) 
    * TR = ((Real,), (Real,)) -> Elliptical Cylinder Mantle (a = r[1][1], b = r[2][1])
    * TR = ((Real, Real),(Real, Real)) -> Elliptical Cone Mantle \n(a_in = r[1][1], a_out = r[1][2], b_in = r[2][1], b_out = r[2][2])
    * Not all are implemented yet

* `φ::TP`: 
    * TP = Nothing <-> Full in φ
    * ...
* `hZ::T`: half hight/length of the cone mantle

* `D`: `:inwards` or `:outwards`: Whethe the normal points inside or outside
"""
@with_kw struct ConeMantle{T,TR,TP,D} <: AbstractCurvedSurfacePrimitive{T}
    r::TR = 1
    φ::TP = nothing
    hZ::T = 1

    # axis::Line{T} = Line{T}(CartesianPoint{T}(zero(T), zero(T), -hZ), CartesianVector{T}(zero(T), zero(T), 2hZ))
    origin::CartesianPoint{T} = zero(CartesianPoint{T})
    rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})
end

radius_at_z(hZ::T, rBot::T, rTop::T, z::T) where {T} = iszero(hZ) ? rBot : rBot + (hZ+z)*(rTop - rBot)/(2hZ) 
radius_at_z(cm::ConeMantle{T,T}, z::T) where {T} = cm.r
radius_at_z(cm::ConeMantle{T,Tuple{T,T}}, z::T) where {T} = radius_at_z(cm.hZ, cm.r[1], cm.r[2], z)

function normal(cm::ConeMantle{T,T}, pt::CartesianPoint{T}) where {T}
    pto = _transform_into_object_coordinate_system(pt, cm)
    cyl = CylindricalPoint(pto)
    return CartesianVector(_transform_into_global_coordinate_system(
            CartesianPoint(CylindricalPoint{T}(cyl.r, cyl.φ, zero(T))), cm))
end
function normal(cm::ConeMantle{T,Tuple{T,T},<:Any,:inwards}, pt::CartesianPoint{T}) where {T}
    pto = _transform_into_object_coordinate_system(pt, cm)
    cyl = CylindricalPoint(pto)
    Δr = cm.r[2] - cm.r[1]
    Δz = 2cm.hZ
    return CartesianVector(_transform_into_global_coordinate_system(
            CartesianPoint(CylindricalPoint{T}(-one(T), cyl.φ, Δr / Δz)), cm))
end
function normal(cm::ConeMantle{T,Tuple{T,T},<:Any,:outwards}, pt::CartesianPoint{T}) where {T}
    pto = _transform_into_object_coordinate_system(pt, cm)
    cyl = CylindricalPoint(pto)
    Δr = cm.r[2] - cm.r[1]
    Δz = 2cm.hZ
    return CartesianVector(_transform_into_global_coordinate_system(
            CartesianPoint(CylindricalPoint{T}( one(T), cyl.φ, -Δr / Δz)), cm))
end

const FullConeMantle{T,D} = ConeMantle{T,Tuple{T,T},Nothing,D} # ugly name but works for now, should just be `ConeMantle`...
const PartialConeMantle{T,D} = ConeMantle{T,Tuple{T,T},Tuple{T,T},D}


extremum(cm::FullConeMantle{T}) where {T} = sqrt(cm.hZ^2 + max(cm.r...)^2)
extremum(cm::PartialConeMantle{T}) where {T} = sqrt(cm.hZ^2 + max(cm.r...)^2)


function lines(sp::FullConeMantle{T}; n = 2) where {T} 
    bot_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T),zero(T), -sp.hZ), sp)
    bot_ellipse = Circle{T}(r = sp.r[1], φ = sp.φ, origin = bot_origin, rotation = sp.rotation)
    top_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T),zero(T), +sp.hZ), sp)
    top_ellipse = Circle{T}(r = sp.r[2], φ = sp.φ, origin = top_origin, rotation = sp.rotation)
    φs = range(T(0), step = T(2π) / n, length = n)
    edges = [ Edge{T}(_transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[1], φ, -sp.hZ)), sp),
                      _transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[2], φ, +sp.hZ)), sp)) for φ in φs ]
    bot_ellipse, top_ellipse, edges
end
function lines(sp::PartialConeMantle{T}; n = 2) where {T} 
    bot_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -sp.hZ), sp)
    bot_ellipse = PartialCircle{T}(r = sp.r[1], φ = sp.φ, origin = bot_origin, rotation = sp.rotation)
    top_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), +sp.hZ), sp)
    top_ellipse = PartialCircle{T}(r = sp.r[2], φ = sp.φ, origin = top_origin, rotation = sp.rotation)
    φs = range(sp.φ[1], stop = sp.φ[2], length = n)
    edges = [ Edge{T}(_transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[1], φ, -sp.hZ)), sp),
                      _transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[2], φ, +sp.hZ)), sp)) for φ in φs ]
    bot_ellipse, top_ellipse, edges
end


"""
    intersection(cm::ConeMantle{T,Tuple{T,T}}, l::Line{T}) where {T}

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
    return _transform_into_global_coordinate_system(ints1, cm), 
           _transform_into_global_coordinate_system(ints2, cm)
end

# """
#     intersection(cm::ConeMantle{T,T}, l::Line{T}) where {T}

# The function will always return 2 CartesianPoint's.
# If the line just touches the mantle, the two points will be the same. 
# If the line does not touch the mantle at all, the two points will have NaN's as there coordinates.
# """
# function intersection(cm::ConeMantle{T,T}, l::Line{T}) where {T}
#     obj_l = _transform_into_object_coordinate_system(l, cm) # direction is not normalized
    
#     L1 = obj_l.origin.x
#     L2 = obj_l.origin.y
#     L3 = obj_l.origin.z
#     D1 = obj_l.direction.x
#     D2 = obj_l.direction.y
#     D3 = obj_l.direction.z

#     f1 = D1^2 + D2^2 
#     λ = inv(f1) # f1 is only 0 if obj_l is parallel to the axis of the cone 
#                 # (here eZ -> D1 = D2 = 0)
#                 # We assume here that this is not the case -> 
#                 # We check this in choosing the sample / evaluating dimensions in `paint!` 
#     hZ = cm.hZ
#     R0 = cm.r

#     term1 = (2D1*L1 + 2D2*L2)^2
#     term2 = L1^2 + L2^2 - R0^2
#     term3 = -D1*L1 - D2*L2 
#     term4 = term1 - 4*f1*term2
#     sq::T = term4 < 0 ? T(NaN) : sqrt(term1 - 4*f1*term2) # if this 
    
#     λ1 = λ * (-sq/2 + term3) 
#     λ2 = λ * (+sq/2 + term3)
    
#     ints1 = obj_l.origin + λ1 * obj_l.direction 
#     ints2 = obj_l.origin + λ2 * obj_l.direction 
#     return _transform_into_global_coordinate_system(ints1, cm), 
#            _transform_into_global_coordinate_system(ints2, cm)
# end


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

# function distance_to_surface(point::AbstractCoordinatePoint{T}, c::ConeMantle{T, <:Any, Nothing, <:Any})::T where {T}
#     pcy = CylindricalPoint(point)
#     distance_to_line(PlanarPoint{T}(pcy.r,pcy.z), LineSegment(c))
# end

# function distance_to_surface(point::AbstractCoordinatePoint{T}, c::ConeMantle{T, <:Any, <:AbstractInterval, <:Any})::T where {T}
#     pcy = CylindricalPoint(point)
#     φMin::T, φMax::T, _ = get_φ_limits(c)
#     if _in_φ(pcy, c.φ)
#         return distance_to_line(PlanarPoint{T}(pcy.r,pcy.z), LineSegment(c))
#     else
#         return distance_to_line(CartesianPoint(point), LineSegment(c, _φNear(pcy.φ, φMin, φMax)))
#     end
# end
