
abstract type AbstractPlanarPoint{T} <: StaticArrays.FieldVector{2, T} end
abstract type AbstractPlanarVector{T} <: StaticArrays.FieldVector{2, T} end

struct PlanarPoint{T} <: AbstractPlanarPoint{T}
    u::T
    v::T
end

"""
    struct CartesianPoint{T} <: AbstractCoordinatePoint{T, Cartesian}

Describes a three-dimensional point in Cartesian coordinates.

## Fields 
* `x`: x-coordinate (in m).
* `y`: y-coordinate (in m).
* `z`: z-coordinate (in m).

See also [`CylindricalPoint`](@ref).
"""
struct CartesianPoint{T} <: AbstractCoordinatePoint{T, Cartesian}
    x::T
    y::T
    z::T
end

zero(PT::Type{<:AbstractCoordinatePoint{T}}) where {T} = PT(zero(T),zero(T),zero(T))

@inline _eq_cyl_r(p::CartesianPoint{T}, r::Real) where {T} = hypot(p.x, p.y) == T(r)

@inline _in_planar_r(p::PlanarPoint, r::Real) = hypot(p.u, p.v) <= r
@inline _in_planar_r(p::PlanarPoint, r::AbstractInterval) = hypot(p.u, p.v) in r

@inline _in_cyl_r(p::CartesianPoint, r::Real) = hypot(p.x, p.y) <= r
@inline _in_cyl_r(p::CartesianPoint, r::AbstractInterval) = hypot(p.x, p.y) in r

@inline _isapprox_φ(p::CartesianPoint{T}, φ::Real) where {T} = isapprox(mod(atan(p.y, p.x), T(2π)), T(φ), atol = geom_atol_zero(T))

@inline _in_φ(p::CartesianPoint{T}, φ::AbstractInterval) where {T} = _in_angular_interval_closed(mod(atan(p.y, p.x), T(2π)), φ)

@inline _in_planar_α(p::PlanarPoint{T}, α::AbstractInterval) where {T} = _in_angular_interval_closed(mod(atan(p.v, p.u), T(2π)), α)

@inline _isapprox_x(p::CartesianPoint{T}, x::Real) where {T} = isapprox(p.x, x, atol = geom_atol_zero(T))

@inline _in_x(p::CartesianPoint, x::Real) = abs(p.x) <= x
@inline _in_x(p::CartesianPoint, x::AbstractInterval) = p.x in x

@inline _in_planar_u(p::PlanarPoint, u::AbstractInterval) = p.u in u

@inline _isapprox_y(p::CartesianPoint{T}, y::Real) where {T} = isapprox(p.y, y, atol = geom_atol_zero(T))

@inline _in_y(p::CartesianPoint, y::Real) = abs(p.y) <= y
@inline _in_y(p::CartesianPoint, y::AbstractInterval) = p.y in y

@inline _in_planar_v(p::PlanarPoint, v::AbstractInterval) = p.v in v

@inline _isapprox_z(p::CartesianPoint{T}, z::Real) where {T} = isapprox(p.z, z, atol = geom_atol_zero(T))

@inline _in_z(p::CartesianPoint, z::Real) = abs(p.z) <= z
@inline _in_z(p::CartesianPoint, z::AbstractInterval) = p.z in z

@inline _isapprox_sph_r(p::CartesianPoint{T}, radius::Real) where {T} = isapprox(hypot(p.x, p.y, p.z), radius, atol = geom_atol_zero(T))

@inline _in_sph_r(p::CartesianPoint, radius::Real) = hypot(p.x, p.y, p.z) <= radius
@inline _in_sph_r(p::CartesianPoint, radius::AbstractInterval) = hypot(p.x, p.y, p.z) in radius

@inline _isapprox_torr_r_tube(p::CartesianPoint{T}, r_torus::Real, r_tube::Real, z_torus::Real) where {T} = isapprox(hypot(hypot(p.x, p.y) - r_torus, p.z - z_torus), T(r_tube), atol = geom_atol_zero(T))

@inline _in_torr_r_tube(p::CartesianPoint, r_torus::Real, r_tube::Real, z_torus::Real) = hypot(hypot(p.x, p.y) - r_torus, p.z - z_torus) <= r_tube
@inline _in_torr_r_tube(p::CartesianPoint, r_torus::Real, r_tube::AbstractInterval, z_torus::Real) = hypot(hypot(p.x, p.y) - r_torus, p.z - z_torus) in r_tube

@inline _in_torr_θ(p::CartesianPoint{T}, r_torus::Real, θ::AbstractInterval, z_torus::Real) where {T} = _in_angular_interval_closed(mod(atan(p.z - z_torus, hypot(p.x, p.y) - r_torus), T(2π)), θ)

"""
    struct CylindricalPoint{T} <: AbstractCoordinatePoint{T, Cylindrical}

Describes a three-dimensional point in cylindrical coordinates. 

## Fields
* `r`: Radius (in m).
* `φ`: Polar angle (in rad).
* `z`: `z`-coordinate (in m).

!!! note 
    `φ == 0` corresponds to the `x`-axis in the Cartesian coordinate system.
    
See also [`CartesianPoint`](@ref).
"""
struct CylindricalPoint{T} <: AbstractCoordinatePoint{T, Cylindrical}
    r::T
    φ::T
    z::T
    CylindricalPoint{T}(r::T, φ::T, z::T) where {T} = new(r, mod(φ,T(2π)), z)
    CylindricalPoint{T}(r::Real, φ::Real, z::Real) where {T} = new(T(r), mod(T(φ),T(2π)), T(z))
end

function CylindricalPoint(pt::CartesianPoint{T})::CylindricalPoint{T} where {T}
    return CylindricalPoint{T}(hypot(pt.x, pt.y), atan(pt.y, pt.x), pt.z)
end

function CartesianPoint(pt::CylindricalPoint{T})::CartesianPoint{T} where {T}
    sφ::T, cφ::T = sincos(pt.φ)
    return CartesianPoint{T}(pt.r * cφ, pt.r * sφ, pt.z)
end

@inline CylindricalPoint(pt::CylindricalPoint) = pt
@inline CartesianPoint(pt::CartesianPoint) = pt

@inline _convert_point(pt::AbstractCoordinatePoint, ::Type{Cylindrical}) = CylindricalPoint(pt)
@inline _convert_point(pt::AbstractCoordinatePoint, ::Type{Cartesian}) = CartesianPoint(pt)

@inline _eq_cyl_r(p::CylindricalPoint{T}, r::Real) where {T} = p.r == T(r)

@inline _in_cyl_r(p::CylindricalPoint, r::Real) = p.r <= r
@inline _in_cyl_r(p::CylindricalPoint, r::AbstractInterval) = p.r in r

@inline _isapprox_φ(p::CylindricalPoint{T}, φ::Real) where {T} = isapprox(p.φ, mod(T(φ), T(2π)), atol = geom_atol_zero(T))

@inline _in_φ(p::CylindricalPoint, φ::AbstractInterval) = _in_angular_interval_closed(p.φ, φ)

@inline _isapprox_x(p::CylindricalPoint{T}, x::Real) where {T} = isapprox(p.r * cos(p.φ), x, atol = geom_atol_zero(T))

@inline _in_x(p::CylindricalPoint, x::Real) = abs(p.r * cos(p.φ)) <= x
@inline _in_x(p::CylindricalPoint, x::AbstractInterval) = p.r * cos(p.φ) in x

@inline _isapprox_y(p::CylindricalPoint{T}, y::Real) where {T} = isapprox(p.r * sin(p.φ), y, atol = geom_atol_zero(T))

@inline _in_y(p::CylindricalPoint, y::Real) = abs(p.r * sin(p.φ)) <= y
@inline _in_y(p::CylindricalPoint, y::AbstractInterval) = p.r * sin(p.φ) in y

@inline _isapprox_z(p::CylindricalPoint{T}, z::Real) where {T} = isapprox(p.z, z, atol = geom_atol_zero(T))

@inline _in_z(p::CylindricalPoint, z::Real) = abs(p.z) <= z
@inline _in_z(p::CylindricalPoint, z::AbstractInterval) = p.z in z

@inline _isapprox_sph_r(p::CylindricalPoint{T}, radius::Real) where {T} = isapprox(hypot(p.r, p.z), radius, atol = geom_atol_zero(T))

@inline _in_sph_r(p::CylindricalPoint, radius::Real) = hypot(p.r, p.z) <= radius
@inline _in_sph_r(p::CylindricalPoint, radius::AbstractInterval) = hypot(p.r, p.z) in radius

@inline _isapprox_torr_r_tube(p::CylindricalPoint{T}, r_torus::Real, r_tube::Real, z_torus::Real) where {T} = isapprox(hypot(p.r - r_torus, p.z - z_torus), T(r_tube), atol = geom_atol_zero(T))
@inline _in_torr_r_tube(p::CylindricalPoint, r_torus::Real, r_tube::Real, z_torus::Real) = hypot(p.r - r_torus, p.z - z_torus) <= r_tube
@inline _in_torr_r_tube(p::CylindricalPoint, r_torus::Real, r_tube::AbstractInterval, z_torus::Real) = hypot(p.r - r_torus, p.z - z_torus) in r_tube

@inline _in_torr_θ(p::CylindricalPoint{T}, r_torus::Real, θ::AbstractInterval, z_torus::Real) where {T} = _in_angular_interval_closed(mod(atan(p.z - z_torus, p.r - r_torus), T(2π)), θ)

function _Δφ(φ1::T, φ2::T)::T where {T}
    δφ = mod(φ2 - φ1, T(2π))
    min(δφ, T(2π) - δφ)
end

_φNear(φ::Real, φMin::T, φMax::T) where {T} = _Δφ(T(φ),φMin) ≤ _Δφ(T(φ),φMax) ? φMin : φMax
