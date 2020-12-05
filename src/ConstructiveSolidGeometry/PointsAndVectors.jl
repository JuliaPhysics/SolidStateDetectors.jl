abstract type AbstractCoordinatePoint{T, S} <: StaticArrays.FieldVector{3, T} end
abstract type AbstractCoordinateVector{T, S} <: StaticArrays.FieldVector{3, T} end


"""
    struct CartesianPoint{T} <: AbstractCoordinatePoint{T, Cartesian}

* `x`: x-coordinate in meter
* `y`: y-coordinate in meter
* `z`: z-coordinate in meter
"""
struct CartesianPoint{T} <: AbstractCoordinatePoint{T, Cartesian}
    x::T
    y::T
    z::T
end

@inline (*)(r::RotMatrix{3,T,TT}, p::CartesianPoint{T}) where {T, TT} = r.mat * p
@inline (*)(r::RotMatrix{3,T,TT}, vp::Vector{CartesianPoint{T}}) where {T, TT} = (r,) .* vp
@inline (*)(r::RotMatrix{3,T,TT}, vvp::Vector{Vector{CartesianPoint{T}}}) where {T, TT} = (r,) .* vvp

@inline (*)(p::CartesianPoint{T}, v::SVector{3,T}) where {T} = p .* v
@inline (*)(vp::Vector{CartesianPoint{T}}, v::SVector{3,T}) where {T} = vp .* (v,)
@inline (*)(vvp::Vector{Vector{CartesianPoint{T}}}, v::SVector{3,T}) where {T} = vvp .* (v,)

@inline _in_cyl_r(p::CartesianPoint, r::Real) = hypot(p.x, p.y) <= r
@inline _in_cyl_r(p::CartesianPoint, r::AbstractInterval) = hypot(p.x, p.y) in r

@inline _in_φ(p::CartesianPoint{T}, φ::AbstractInterval) where {T} = mod(atan(p.y, p.x), T(2π)) in φ    

@inline _in_x(p::CartesianPoint, x::Real) = abs(p.x) <= x
@inline _in_x(p::CartesianPoint, x::AbstractInterval) = p.x in x

@inline _in_y(p::CartesianPoint, y::Real) = abs(p.y) <= y
@inline _in_y(p::CartesianPoint, y::AbstractInterval) = p.y in y

@inline _in_z(p::CartesianPoint, z::Real) = abs(p.z) <= z
@inline _in_z(p::CartesianPoint, z::AbstractInterval) = p.z in z

@inline _in_sph_r(p::CartesianPoint, radius::Real) = hypot(p.x, p.y, p.z) <= radius
@inline _in_sph_r(p::CartesianPoint, radius::AbstractInterval) = hypot(p.x, p.y, p.z) in radius

"""
    struct CylindricalPoint{T} <: AbstractCoordinatePoint{T, Cylindrical}

* `r`: Radius in meter
* `φ`: Polar angle in radians. φ == 0 <=> Parallel to x-axis of cartesian coordinate system."
* `z`: z-coordinate in meter
"""
struct CylindricalPoint{T} <: AbstractCoordinatePoint{T, Cylindrical}
    r::T
    φ::T
    z::T
    CylindricalPoint{T}(r::T, φ::T, z::T) where {T} = new(r, mod(φ,T(2π)), z)
end

function CylindricalPoint(pt::CartesianPoint{T})::CylindricalPoint{T} where {T}
    return CylindricalPoint{T}(hypot(pt.x, pt.y), atan(pt.y, pt.x), pt.z)
end

function CartesianPoint(pt::CylindricalPoint{T})::CartesianPoint{T} where {T}
    sφ::T, cφ::T = sincos(pt.φ)
    return CartesianPoint{T}(pt.r * cφ, pt.r * sφ, pt.z)
end

@inline _in_cyl_r(p::CylindricalPoint, r::Real) = p.r <= r
@inline _in_cyl_r(p::CylindricalPoint, r::AbstractInterval) = p.r in r

@inline _in_φ(p::CylindricalPoint, φ::AbstractInterval) = p.φ in φ

@inline _in_x(p::CylindricalPoint, x::Real) = abs(p.r * cos(p.φ)) <= x
@inline _in_x(p::CylindricalPoint, x::AbstractInterval) = p.r * cos(p.φ) in x

@inline _in_y(p::CylindricalPoint, y::Real) = abs(p.r * sin(p.φ)) <= y
@inline _in_y(p::CylindricalPoint, y::AbstractInterval) = p.r * sin(p.φ) in y

@inline _in_z(p::CylindricalPoint, z::Real) = abs(p.z) <= z
@inline _in_z(p::CylindricalPoint, z::AbstractInterval) = p.z in z

@inline _in_sph_r(p::CylindricalPoint, radius::Real) = hypot(p.r, p.z) <= radius
@inline _in_sph_r(p::CylindricalPoint, radius::AbstractInterval) = hypot(p.r, p.z) in radius


"""
    struct CartesianVector{T} <: AbstractCoordinateVector{T, Cartesian}

* `x`: x-component in meter
* `y`: y-component in meter
* `z`: z-component in meter
"""
struct CartesianVector{T} <: AbstractCoordinateVector{T, Cartesian}
    x::T
    y::T
    z::T
end

@inline (+)(vp::Vector{CartesianPoint{T}}, v::CartesianVector{T}) where {T} = vp .+ (v,)
@inline (+)(vvp::Vector{Vector{CartesianPoint{T}}}, v::CartesianVector{T}) where {T} = vvp .+ (v,)

"""
    struct CylindricalVector{T} <: AbstractCoordinateVector{T, Cylindrical}

* `r`: Radius in meter
* `φ`: Polar angle in radians. φ == 0 <=> Parallel to x-axis of cartesian coordinate system."
* `z`: z-coordinate in meter
"""
struct CylindricalVector{T} <: AbstractCoordinateVector{T, Cylindrical}
    r::T
    φ::T
    z::T
end