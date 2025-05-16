# abstract type AbstractPlanarPoint{T} <: StaticArrays.FieldVector{2, T} end
# abstract type AbstractPlanarVector{T} <: StaticArrays.FieldVector{2, T} end
# 
# struct PlanarPoint{T} <: AbstractPlanarPoint{T}
#     u::T
#     v::T
# end

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

#Type promotion happens here
function CartesianPoint(x::TX, y::TY, z::TZ) where {TX<:Real,TY<:Real,TZ<:Real}
    eltypes = _csg_get_promoted_eltype.((TX,TY,TZ))
    T = float(promote_type(eltypes...))
    CartesianPoint{T}(T(x),T(y),T(z))
end

function CartesianPoint(;
    x = 0,
    y = 0,
    z = 0
)
    CartesianPoint(x,y,z)
end

function CartesianPoint{T}(;
    x = 0,
    y = 0,
    z = 0
) where {T}
    CartesianPoint{T}(T(x),T(y),T(z))
end

function Base.convert(::Type{CartesianPoint{T}}, pt::CartesianPoint{U}) where {T,U}
    return CartesianPoint{T}(convert(T, pt.x), convert(T, pt.y), convert(T, pt.z))
end


Base.:(==)(a::CartesianPoint, b::CartesianPoint) = a.x == b.x && a.y == b.y && a.z == b.z

function Base.isapprox(a::CartesianPoint, b::CartesianPoint; kwargs...)
    return isapprox(a.x, b.x; kwargs...) && isapprox(a.y, b.y; kwargs...) && isapprox(a.z, b.z; kwargs...)
end


@inline Base.:(+)(pt::CartesianPoint, v::CartesianVector) = CartesianPoint(pt.x + v.x, pt.y + v.y, pt.z + v.z)
@inline Base.:(-)(pt::CartesianPoint, v::CartesianVector) = CartesianPoint(pt.x - v.x, pt.y - v.y, pt.z - v.z)
@inline Base.:(-)(a::CartesianPoint, b::CartesianPoint) = CartesianVector(a.x - b.x, a.y - b.y, a.z - b.z)


Base.:(*)(A::StaticMatrix{3,3}, pt::CartesianPoint) = _ascartpoint(A * _asvector(pt))

@inline _asvector(pt::CartesianPoint{T}) where {T} = SVector(pt.x, pt.y, pt.z)
@inline _ascartpoint(v::StaticVector{3}) = CartesianPoint(v[1], v[2], v[3])



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
    CylindricalPoint{T}(r::T, φ::T, z::T) where {T<:AbstractFloat} = new(r, mod(φ,T(2π)), z)
    CylindricalPoint{T}(r::Real, φ::Real, z::Real) where {T} = new(T(r), mod(T(φ),T(2π)), T(z))
end

function CylindricalPoint(r::TR, φ::TP, z::TZ) where {TR<:Real,TP<:Real,TZ<:Real}
    eltypes = _csg_get_promoted_eltype.((TR,TP,TZ))
    T = float(promote_type(eltypes...))
    CylindricalPoint{T}(T(r),T(φ),T(z))
end

function CylindricalPoint(;
    r = 0,
    φ = 0,
    z = 0
)
    CylindricalPoint(r,φ,z)
end

function CylindricalPoint{T}(;
    r = 0,
    φ = 0,
    z = 0
) where {T}
    CylindricalPoint{T}(T(r),T(φ),T(z))
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


function Base.convert(::Type{CylindricalPoint{T}}, pt::CylindricalPoint{U}) where {T,U}
    return CylindricalPoint{T}(convert(T, pt.r), convert(T, pt.φ), convert(T, pt.z))
end


Base.:(==)(a::CylindricalPoint, b::CylindricalPoint) = a.r == b.r && a.φ == b.φ && a.z == b.z

function Base.isapprox(a::CylindricalPoint, b::CylindricalPoint; kwargs...)
    return isapprox(a.r, b.r; kwargs...) && isapprox(a.φ, b.φ; kwargs...) && isapprox(a.z, b.z; kwargs...)
end


@inline Base.:(+)(pt::CylindricalPoint, v::CylindricalVector) = CylindricalPoint(pt.r + v.r, pt.φ + v.φ, pt.z + v.z)

@inline Base.:(-)(pt::CylindricalPoint, v::CylindricalVector) = CylindricalPoint(pt.r - v.r, pt.φ - v.φ, pt.z - v.z)

@inline Base.:(-)(a::CylindricalPoint, b::CylindricalPoint) = CylindricalVector(a.r - b.r, a.φ - b.φ, a.z - b.z)


# function _Δφ(φ1::T, φ2::T)::T where {T}
#     δφ = mod(φ2 - φ1, T(2π))
#     min(δφ, T(2π) - δφ)
# end
# 
# _φNear(φ::Real, φMin::T, φMax::T) where {T} = _Δφ(T(φ),φMin) ≤ _Δφ(T(φ),φMax) ? φMin : φMax
