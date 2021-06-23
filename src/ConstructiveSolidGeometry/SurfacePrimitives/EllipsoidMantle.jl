@with_kw struct EllipsoidMantle{T,RT,PT,TT} <: AbstractCurvedSurfacePrimitive{T}
    r::RT = 1
    φ::PT = nothing
    θ::TT = nothing

    origin::CartesianPoint{T} = zero(CartesianPoint{T})
    rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})
end

const FullSphereMantle{T} = EllipsoidMantle{T,T,Nothing,Nothing}
const FullEllipsoidMantle{T} = EllipsoidMantle{T,NTuple{3,T},Nothing,Nothing}

# get_r_limits(s::SphereMantle{T}) where {T} = (_left_radial_interval(s.r), _right_radial_interval(s.r))
# get_φ_limits(s::SphereMantle{T}) where {T} = (T(0), T(2π), true)
# get_z_limits(s::SphereMantle{T}) where {T} = (-_right_linear_interval(s.r), _right_linear_interval(s.r))

# in(p::AbstractCoordinatePoint, s::SphereMantle) = _isapprox_sph_r(p, s.r)

# #=
# function sample(s::SphereMantle{T}, step::Real)::Vector{CylindricalPoint{T}} where {T}
#     samples = [
#         CylindricalPoint{T}(sqrt(s.r^2 - z^2),φ,z)
#         for z in in -s.r:step:s.r
#         for φ in 0:step/s.r:2π
#     ]
# end
# =#

# function sample(s::SphereMantle{T}, Nsamps::NTuple{3,Int})::Vector{CylindricalPoint{T}} where {T}
#     samples = [
#         CylindricalPoint{T}(sqrt(s.r^2 - z^2),φ,z)
#         for z in (Nsamps[3] ≤ 1 ? s.r : range(-s.r, s.r, length = Nsamps[3]))
#         for φ in (Nsamps[2] ≤ 1 ? 0 : range(0, 2π, length = Nsamps[2]))
#     ]
# end

# function sample(s::SphereMantle{T}, g::CylindricalTicksTuple{T})::Vector{CylindricalPoint{T}} where {T}
#     samples = [
#         CylindricalPoint{T}((s.r^2 - z^2 < 0 ? 0 : sqrt(s.r^2 - z^2)),φ,z)
#         for z in get_z_ticks(s, g)
#         for φ in get_φ_ticks(s, g)
#     ]
# end

# function _get_x_at_z(s::SphereMantle{T}, g::CartesianTicksTuple{T}, z::T) where {T}
#     if s.r < abs(z) return [0] end
#     R::T = sqrt(s.r^2 - z^2)
#     x_from_y::Vector{T} = sqrt.(R^2 .- filter(y -> abs(y) <= R, g.y).^2)
#     _get_ticks(sort!(vcat(g.x, x_from_y, -x_from_y)), -R, R)
# end

# function _get_y_at_z(s::SphereMantle{T}, z::T, x::T) where {T}
#     tmp::T = s.r^2 - hypot(x,z)^2
#     if tmp < 0 return (0,) end
#     (-sqrt(tmp), sqrt(tmp))
# end

# function sample(s::SphereMantle{T}, g::CartesianTicksTuple{T})::Vector{CartesianPoint{T}} where {T}
#     samples = [
#         CartesianPoint{T}(x,y,z)
#         for z in get_z_ticks(s, g)
#         for x in _get_x_at_z(s, g, z)
#         for y in _get_y_at_z(s, z, x)
#     ]
# end

# distance_to_surface(point::CylindricalPoint{T}, s::SphereMantle{T}) where {T} = abs(hypot(point.r, point.z) - s.r)

# distance_to_surface(point::CartesianPoint{T}, s::SphereMantle{T}) where {T} = abs(norm(point) - s.r)
