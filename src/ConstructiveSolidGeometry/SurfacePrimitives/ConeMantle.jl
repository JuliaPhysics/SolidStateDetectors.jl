"""
    struct ConeMantle{T,RT,TP} <: AbstractSurfacePrimitive{T}

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
* `zH::T`: half hight/length of the cone mantle

* `axis::Line{T}`: The axis of the Cone mantle. Going from the bottom center point to the top center point. 
"""
@with_kw struct ConeMantle{T,RT,TP} <: AbstractCurvedSurfacePrimitive{T}
    r::RT = 1
    φ::TP = nothing
    hZ::T = 1 # maybe we don't need this. I will leave it for now...

    # axis::Line{T} = Line{T}(CartesianPoint{T}(zero(T), zero(T), -hZ), CartesianVector{T}(zero(T), zero(T), 2hZ))
    origin::CartesianPoint{T} = zero(CartesianPoint{T})
    rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})
end


const CylinderMantle{T} = ConeMantle{T,T,Nothing}
const PartialCylinderMantle{T} = ConeMantle{T,T,Tuple{T,T}}

const FullConeMantle{T} = ConeMantle{T,Tuple{T,T},Nothing} # ugly name but works for now, should just be `ConeMantle`...
const PartialConeMantle{T} = ConeMantle{T,Tuple{T,T},Tuple{T,T}}

extremum(cm::CylinderMantle{T}) where {T} = sqrt(cm.hZ^2 + cm.r^2)
extremum(cm::PartialCylinderMantle{T}) where {T} = sqrt(cm.hZ^2 + cm.r^2)

extremum(cm::FullConeMantle{T}) where {T} = sqrt(cm.hZ^2 + max(es.r[1], es.r[2])^2)
extremum(cm::PartialConeMantle{T}) where {T} = sqrt(cm.hZ^2 + max(es.r[1], es.r[2])^2)

function lines(sp::CylinderMantle{T}) where {T} 
    bot_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T),  sp.hZ), sp)
    bot_ellipse = Circle{T}(r = sp.r, φ = sp.φ, origin = bot_origin, rotation = sp.rotation)
    top_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -sp.hZ), sp)
    top_ellipse = Circle{T}(r = sp.r, φ = sp.φ, origin = top_origin, rotation = sp.rotation)
    bot_ellipse, top_ellipse
end
function lines(sp::PartialCylinderMantle{T}) where {T} 
    bot_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T),  sp.hZ), sp)
    bot_ellipse = PartialCircle{T}(r = sp.r, φ = sp.φ, origin = bot_origin, rotation = sp.rotation)
    top_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -sp.hZ), sp)
    top_ellipse = PartialCircle{T}(r = sp.r, φ = sp.φ, origin = top_origin, rotation = sp.rotation)
    p_bot_l = _transform_into_global_coordinate_system(CartesianPoint{T}(sp.r, sp.φ[1], -sp.hZ), sp)
    p_bot_r = _transform_into_global_coordinate_system(CartesianPoint{T}(sp.r, sp.φ[2], -sp.hZ), sp)
    p_top_l = _transform_into_global_coordinate_system(CartesianPoint{T}(sp.r, sp.φ[1],  sp.hZ), sp)
    p_top_r = _transform_into_global_coordinate_system(CartesianPoint{T}(sp.r, sp.φ[2],  sp.hZ), sp)
    edge_l = Edge{T}(p_bot_l, p_top_l)
    edge_r = Edge{T}(p_bot_r, p_top_r)
    bot_ellipse, top_ellipse, edge_l, edge_r
end
function lines(sp::FullConeMantle{T}) where {T} 
    bot_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T),zero(T), sp.hZ), sp)
    bot_ellipse = Circle{T}(r = sp.r[1], φ = sp.φ, origin = bot_origin, rotation = sp.rotation)
    top_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T),zero(T),-sp.hZ), sp)
    top_ellipse = Circle{T}(r = sp.r[2], φ = sp.φ, origin = top_origin, rotation = sp.rotation)
    bot_ellipse, top_ellipse
end
function lines(sp::PartialConeMantle{T}) where {T} 
    bot_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T),  sp.hZ), sp)
    bot_ellipse = PartialCircle{T}(r = sp.r[1], φ = sp.φ, origin = bot_origin, rotation = sp.rotation)
    top_origin = _transform_into_global_coordinate_system(CartesianPoint{T}(zero(T), zero(T), -sp.hZ), sp)
    top_ellipse = PartialCircle{T}(r = sp.r[2], φ = sp.φ, origin = top_origin, rotation = sp.rotation)
    p_bot_l = _transform_into_global_coordinate_system(CartesianPoint{T}(sp.r[1], sp.φ[1], -sp.hZ), sp)
    p_bot_r = _transform_into_global_coordinate_system(CartesianPoint{T}(sp.r[1], sp.φ[2], -sp.hZ), sp)
    p_top_l = _transform_into_global_coordinate_system(CartesianPoint{T}(sp.r[2], sp.φ[1],  sp.hZ), sp)
    p_top_r = _transform_into_global_coordinate_system(CartesianPoint{T}(sp.r[2], sp.φ[2],  sp.hZ), sp)
    edge_l = Edge{T}(p_bot_l, p_top_l)
    edge_r = Edge{T}(p_bot_r, p_top_r)
    bot_ellipse, top_ellipse, edge_l, edge_r
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
        sq::T = term4 < 0 ? T(NaN) : sqrt(term4) 
    
        λ1 = λ * (-sq/2 + term3) 
        λ2 = λ * (+sq/2 + term3)
        λ1, λ2    
    end
    ints1 = obj_l.origin + λ1 * obj_l.direction 
    ints2 = obj_l.origin + λ2 * obj_l.direction 
    return _transform_into_global_coordinate_system(ints1, cm), 
           _transform_into_global_coordinate_system(ints2, cm)
end

"""
    intersection(cm::CylinderMantle{T}, l::Line{T}) where {T}

The function will always return 2 CartesianPoint's.
If the line just touches the mantle, the two points will be the same. 
If the line does not touch the mantle at all, the two points will have NaN's as there coordinates.
"""
function intersection(cm::CylinderMantle{T}, l::Line{T}) where {T}
    obj_l = _transform_into_object_coordinate_system(l, cm) # direction is not normalized
    
    L1 = obj_l.origin.x
    L2 = obj_l.origin.y
    L3 = obj_l.origin.z
    D1 = obj_l.direction.x
    D2 = obj_l.direction.y
    D3 = obj_l.direction.z

    f1 = D1^2 + D2^2 
    λ = inv(f1) # f1 is only 0 if obj_l is parallel to the axis of the cone 
                # (here eZ -> D1 = D2 = 0)
                # We assume here that this is not the case -> 
                # We check this in choosing the sample / evaluating dimensions in `paint!` 
    hZ = cm.hZ
    R0 = cm.r

    term1 = (2D1*L1 + 2D2*L2)^2
    term2 = L1^2 + L2^2 - R0^2
    term3 = -D1*L1 - D2*L2 
    term4 = term1 - 4*f1*term2
    sq::T = term4 < 0 ? T(NaN) : sqrt(term1 - 4*f1*term2) # if this 
    
    λ1 = λ * (-sq/2 + term3) 
    λ2 = λ * (+sq/2 + term3)
    
    ints1 = obj_l.origin + λ1 * obj_l.direction 
    ints2 = obj_l.origin + λ2 * obj_l.direction 
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


# function ConeMantle(c::Cone{T}; rbot = 1, rtop = 1) where {T}
#     r = rbot == rtop ? T(rbot) : (T(rbot), T(rtop))
#     ConeMantle( T, r, c.φ, c.z)
# end

# function ConeMantle(t::Torus{T}; θ = π/2) where {T}
#     r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
#     θ = T(mod(θ,2π))
#     sθ, cθ = sincos(θ)
#     if θ > T(0) && θ < T(π)
#         rbot = geom_round(t.r_torus + r_tubeMin*cθ)
#         rtop = geom_round(t.r_torus + r_tubeMax*cθ)
#         zMin = geom_round(t.z + r_tubeMin*sθ)
#         zMax = geom_round(t.z + r_tubeMax*sθ)
#     elseif θ > T(π) && θ < T(2π)
#         rtop = geom_round(t.r_torus + r_tubeMin*cθ)
#         rbot = geom_round(t.r_torus + r_tubeMax*cθ)
#         zMax = geom_round(t.z + r_tubeMin*sθ)
#         zMin = geom_round(t.z + r_tubeMax*sθ)
#     else
#         @error "Cone Mantle not defined for torroidal cordinate θ = 0 or θ = π. Use Annulus"
#     end
#     r = rbot == rtop ? T(rbot) : (T(rbot), T(rtop))
#     z = zMax == -zMin ? T(zMax) : T(zMin)..T(zMax)
#     ConeMantle( T, r, t.φ, z)
# end

# function ConeMantle(;rbot = 1, rtop = 0, φMin = 0, φMax = 2π, zMin = -1/2, zMax = 1/2)
#     T = float(promote_type(typeof.((rbot, rtop, φMin, φMax, zMin, zMax))...))
#     r = rbot == rtop ? T(rbot) : (T(rbot), T(rtop))
#     φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
#     z = zMax == -zMin ? T(zMax) : T(zMin)..T(zMax)
#     ConeMantle(T, r, φ, z)
# end
# ConeMantle(rbot, rtop, φMin, φMax, zMin, zMax) = ConeMantle(;rbot = rbot, rtop = rtop, φMin = φMin, φMax = φMax, zMin = zMin, zMax = zMax)

# function ConeMantle(rbot::R1, rtop::R2, height::H) where {R1<:Real, R2<:Real, H<:Real}
#     T = float(promote_type(R1, R2, H))
#     ConeMantle( T, (T(rbot), T(rtop)), nothing, T(height)/2)
# end


# get_r_at_z(c::ConeMantle{T, T}, z::Real) where {T} = c.r
# get_r_at_z(c::ConeMantle{T, Tuple{T,T}}, z::Real) where {T} = _get_r_at_z(c.r[1], c.r[2], c.z, z)

# get_r_limits(c::ConeMantle{T, T, <:Any, <:Any}) where {T} = (T(c.r), T(c.r))
# get_r_limits(c::ConeMantle{T, <:Tuple, <:Any, <:Any}) where {T} = c.r

# get_φ_limits(c::ConeMantle{T, <:Any, Nothing, <:Any}) where {T} = (T(0), T(2π), true)
# get_φ_limits(c::ConeMantle{T, <:Any, <:AbstractInterval, <:Any}) where {T} = (c.φ.left, c.φ.right, false)

# get_z_limits(c::ConeMantle{T}) where {T} = (_left_linear_interval(c.z), _right_linear_interval(c.z))

# in(p::AbstractCoordinatePoint, c::ConeMantle{<:Any, <:Any, Nothing, <:Any}) =
#     _in_z(p, c.z) && _eq_cyl_r(p, get_r_at_z(c, p.z))

# in(p::AbstractCoordinatePoint, c::ConeMantle{<:Any, <:Any, <:AbstractInterval, <:Any}) =
#     _in_z(p, c.z) && _in_φ(p, c.φ) && _eq_cyl_r(p, get_r_at_z(c, p.z))

# #=
# function sample(c::ConeMantle{T}, step::Real)::Vector{CylindricalPoint{T}} where {T}
#     φMin::T, φMax::T, _ = get_φ_limits(c)
#     zMin::T, zMax::T = get_z_limits(c)
#     samples = [
#         CylindricalPoint{T}(get_r_at_z(c, z),φ,z)
#         for z in zMin:step:zMax
#         for φ in (get_r_at_z(c, z) == 0 ? φMin : φMin:step/get_r_at_z(c, z):φMax)
#     ]
# end
# =#

# function sample(c::ConeMantle{T}, Nsamps::NTuple{3,Int})::Vector{CylindricalPoint{T}} where {T}
#     φMin::T, φMax::T, _ = get_φ_limits(c)
#     zMin::T, zMax::T = get_z_limits(c)
#     samples = [
#         CylindricalPoint{T}(get_r_at_z(c, z),φ,z)
#         for z in (Nsamps[3] ≤ 1 ? zMin : range(zMin, zMax, length = Nsamps[3]))
#         for φ in (Nsamps[2] ≤ 1 ? φMin : range(φMin, φMax, length = Nsamps[2]))
#     ]
# end

# function sample(c::ConeMantle{T}, g::CylindricalTicksTuple{T})::Vector{CylindricalPoint{T}} where {T}
#     samples = [
#         CylindricalPoint{T}(get_r_at_z(c, z),φ,z)
#         for z in get_z_ticks(c, g)
#         for φ in get_φ_ticks(c, g)
#     ]
# end

# function _get_x_at_z(c::ConeMantle{T}, g::CartesianTicksTuple, z::T) where {T}
#     R::T = get_r_at_z(c, z)
#     x_from_y::Vector{T} = sqrt.(R^2 .- filter(y -> abs(y) <= R, g.y).^2)
#     _get_ticks(sort!(vcat(g.x, x_from_y, -x_from_y)), -R, R)
# end

# function _get_y_at_z(c::ConeMantle{T}, x::T, z::T) where {T}
#     R::T = get_r_at_z(c, z)
#     (-sqrt(R^2-x^2),sqrt(R^2-x^2))
# end

# function sample(c::ConeMantle{T}, g::CartesianTicksTuple{T})::Vector{CartesianPoint{T}} where {T}
#     samples = [
#         CartesianPoint{T}(x,y,z)
#         for z in _get_ticks(g.z, _left_linear_interval(c.z), _right_linear_interval(c.z))
#         for x in _get_x_at_z(c, g, z)
#         for y in _get_y_at_z(c, x, z)
#         if c.φ === nothing || mod(atan(y, x), T(2π)) in c.φ
#     ]
# end

# function LineSegment(c::ConeMantle{T}) where {T}
#     rbot::T, rtop::T = get_r_limits(c)
#     zMin::T, zMax::T = get_z_limits(c)
#     LineSegment(T, PlanarPoint{T}(rbot,zMin), PlanarPoint{T}(rtop,zMax))
# end

# function LineSegment(c::ConeMantle{T}, φ::Real) where {T}
#     rbot::T, rtop::T = get_r_limits(c)
#     zMin::T, zMax::T = get_z_limits(c)
#     sφ::T, cφ::T = sincos(φ)
#     LineSegment(T, CartesianPoint{T}(rbot*cφ,rbot*sφ,zMin), CartesianPoint{T}(rtop*cφ,rtop*sφ,zMax))
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
