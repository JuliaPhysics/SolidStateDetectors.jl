"""
    struct ConeMantle{T,RT,TP} <: AbstractSurfacePrimitive{T}

T: Type of values, e.g. Float64

* `r::TR`: 
    * TR = Real -> Tube Mantle (a = b = r)
    * TR = (Real, Real) -> Cone Mantle (r_bot = r[1], r_top = r[2])
    * TR = ((Real,), (Real,)) -> Elliptical Tube Mantle (a = r[1][1], b = r[2][1])
    * TR = ((Real, Real),(Real, Real)) -> Elliptical Cone Mantle \n(a_in = r[1][1], a_out = r[1][2], b_in = r[2][1], b_out = r[2][2])
    * Not all are implemented yet

* `φ::TP`: 
    * TP = Nothing <-> Full in φ
    * ...
* `zH::T`: half hight/length of the cone mantle

* `axis::Line{T}`: The axis of the Cone mantle. Going from the bottom center point to the top center point. 
"""
@with_kw struct ConeMantle{T,RT,TP} <: AbstractSurfacePrimitive{T}
    r::RT = 1
    φ::TP = nothing
    hZ::T = 1 # maybe we don't need this. I will leave it for now...

    # axis::Line{T} = Line{T}(CartesianPoint{T}(zero(T), zero(T), -hZ), CartesianVector{T}(zero(T), zero(T), 2hZ))
    origin::CartesianPoint{T} = zero(CartesianPoint{T})
    rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})
end


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
