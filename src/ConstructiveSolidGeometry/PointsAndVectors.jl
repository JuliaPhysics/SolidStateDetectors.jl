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

@inline _eq_cyl_r(p::CartesianPoint{T}, r::Real) where {T} = hypot(p.x, p.y) == T(r)

@inline _in_cyl_r(p::CartesianPoint, r::Real) = hypot(p.x, p.y) <= r
@inline _in_cyl_r(p::CartesianPoint, r::AbstractInterval) = hypot(p.x, p.y) in r

@inline _eq_φ(p::CartesianPoint{T}, φ::Real) where {T} = mod(atan(p.y, p.x), T(2π)) == T(φ)

@inline _in_φ(p::CartesianPoint{T}, φ::AbstractInterval) where {T} = mod(atan(p.y, p.x), T(2π)) in φ

@inline _in_x(p::CartesianPoint, x::Real) = abs(p.x) <= x
@inline _in_x(p::CartesianPoint, x::AbstractInterval) = p.x in x

@inline _in_y(p::CartesianPoint, y::Real) = abs(p.y) <= y
@inline _in_y(p::CartesianPoint, y::AbstractInterval) = p.y in y

@inline _eq_z(p::CartesianPoint{T}, z::Real) where {T} = p.z == T(z)

@inline _in_z(p::CartesianPoint, z::Real) = abs(p.z) <= z
@inline _in_z(p::CartesianPoint, z::AbstractInterval) = p.z in z

@inline _in_sph_r(p::CartesianPoint, radius::Real) = hypot(p.x, p.y, p.z) <= radius
@inline _in_sph_r(p::CartesianPoint, radius::AbstractInterval) = hypot(p.x, p.y, p.z) in radius

@inline _eq_torr_r_tube(p::CartesianPoint{T}, r_torus::Real, r_tube::Real) where {T} = T(hypot(hypot(p.x, p.y) - r_torus, p.z)) == T(r_tube)

@inline _in_torr_r_tube(p::CartesianPoint, r_torus::Real, r_tube::Real) = hypot(hypot(p.x, p.y) - r_torus, p.z) <= r_tube
@inline _in_torr_r_tube(p::CartesianPoint, r_torus::Real, r_tube::AbstractInterval) = hypot(hypot(p.x, p.y) - r_torus, p.z) in r_tube

@inline _in_torr_θ(p::CartesianPoint{T}, r_torus::Real, θ::AbstractInterval) where {T} = mod(atan(p.z, hypot(p.x, p.y) - r_torus), T(2π)) in θ

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

@inline _eq_φ(p::CylindricalPoint{T}, φ::Real) where {T} = p.φ == mod(T(φ), T(2π))

@inline _in_φ(p::CylindricalPoint, φ::AbstractInterval) = p.φ in φ

@inline _in_x(p::CylindricalPoint, x::Real) = abs(p.r * cos(p.φ)) <= x
@inline _in_x(p::CylindricalPoint, x::AbstractInterval) = p.r * cos(p.φ) in x

@inline _in_y(p::CylindricalPoint, y::Real) = abs(p.r * sin(p.φ)) <= y
@inline _in_y(p::CylindricalPoint, y::AbstractInterval) = p.r * sin(p.φ) in y

@inline _eq_z(p::CylindricalPoint{T}, z::Real) where {T} = p.z == T(z)

@inline _in_z(p::CylindricalPoint, z::Real) = abs(p.z) <= z
@inline _in_z(p::CylindricalPoint, z::AbstractInterval) = p.z in z

@inline _in_sph_r(p::CylindricalPoint, radius::Real) = hypot(p.r, p.z) <= radius
@inline _in_sph_r(p::CylindricalPoint, radius::AbstractInterval) = hypot(p.r, p.z) in radius

@inline _eq_torr_r_tube(p::CylindricalPoint{T}, r_torus::Real, r_tube::Real) where {T} = T(hypot(p.r - r_torus, p.z)) == T(r_tube)
@inline _in_torr_r_tube(p::CylindricalPoint, r_torus::Real, r_tube::Real) = hypot(p.r - r_torus, p.z) <= r_tube
@inline _in_torr_r_tube(p::CylindricalPoint, r_torus::Real, r_tube::AbstractInterval) = hypot(p.r - r_torus, p.z) in r_tube

@inline _in_torr_θ(p::CylindricalPoint{T}, r_torus::Real, θ::AbstractInterval) where {T} = mod(atan(p.z, p.r - r_torus), T(2π)) in θ

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


@inline rotate(p::CartesianPoint{T}, r::RotMatrix{3,T,TT}) where {T, TT} = r.mat * p
@inline rotate(p::CylindricalPoint{T}, r::RotMatrix{3,T,TT}) where {T, TT} = CylindricalPoint(rotate(CartesianPoint(p), r))
@inline rotate!(vp::Vector{<:AbstractCoordinatePoint{T}}, r::RotMatrix{3,T,TT}) where {T, TT} = begin for i in eachindex(vp) vp[i] = rotate(vp[i], r) end; vp end
@inline rotate!(vvp::Vector{<:Vector{<:AbstractCoordinatePoint{T}}}, r::RotMatrix{3,T,TT}) where {T, TT} = begin for i in eachindex(vvp) rotate!(vvp[i], r) end; vvp end

@inline scale(p::CartesianPoint{T}, v::SVector{3,T}) where {T} = p .* v
@inline scale(p::CylindricalPoint{T}, v::SVector{3,T}) where {T} = CylindricalPoint(scale(CartesianPoint(p),v))
@inline scale!(vp::Vector{<:AbstractCoordinatePoint{T}}, v::SVector{3,T}) where {T} = begin for i in eachindex(vp) vp[i] = scale(vp[i], v) end; vp end
@inline scale!(vvp::Vector{<:Vector{<:AbstractCoordinatePoint{T}}}, v::SVector{3,T}) where {T} = begin for i in eachindex(vvp) scale!(vvp[i], v) end; vvp end

@inline translate(p::CartesianPoint{T}, v::CartesianVector{T}) where {T} = p + v
@inline translate(p::CylindricalPoint{T}, v::CartesianVector{T}) where {T} = CylindricalPoint(translate(CartesianPoint(p), v))
@inline translate!(vp::Vector{<:AbstractCoordinatePoint{T}}, v::CartesianVector{T}) where {T} = begin for i in eachindex(vp) vp[i] = translate(vp[i], v) end; vp end
@inline translate!(vvp::Vector{<:Vector{<:AbstractCoordinatePoint{T}}}, v::CartesianVector{T}) where {T} =  begin for i in eachindex(vvp) translate!(vvp[i], v) end; vvp end


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
