"""
    const AbstractCartesianPoint{T} = AbstractCoordinatePoint{T,Cartesian}

Supertype for cartesian point types.
"""
const AbstractCartesianPoint{T} = AbstractCoordinatePoint{T,Cartesian}

Base.:(==)(a::AbstractCartesianPoint, b::AbstractCartesianPoint) = get_x(a) == get_x(b) && get_y(a) == get_y(b) && get_z(a) == get_z(b)

function Base.isapprox(a::AbstractCartesianPoint, b::AbstractCartesianPoint; kwargs...)
    return isapprox(get_x(a), get_x(b); kwargs...) &&
           isapprox(get_y(a), get_y(b); kwargs...) &&
           isapprox(get_z(a), get_z(b); kwargs...)
end

function to_internal_units(x::AbstractCartesianPoint)
    CartesianPoint(to_internal_units(get_x(x)), to_internal_units(get_y(x)), to_internal_units(get_z(x)))
end



# Unitful uses the `*` and `/` operators to combine values with units. But mathematically that's really not an algebraic
# product (not defined for affine points), but a cartesian product, so supporting this should be fine:
Base.:(*)(pt::AbstractCartesianPoint{<:Real}, u::Unitful.Units{<:Any,Unitful.𝐋}) = CartesianPoint(get_x(pt) * u, get_y(pt) * u, get_z(pt) * u)
Base.:(/)(pt::AbstractCartesianPoint{<:Quantity{<:Real, Unitful.𝐋}}, u::Unitful.Units{<:Any,Unitful.𝐋}) = CartesianPoint(get_x(pt) / u, get_y(pt) / u, get_z(pt) / u)


# In that spirit, and to enable constructs like `ustrip(u"mm", pt)` and `NoUnits(pt / u"mm")`, we'll support uconvert/ustrip
# for CartesianPoint. Unitful doesn't encourage defining those for collections, but we'll view cartesian points as single
# mathematical objects in regard to units:
function Unitful.uconvert(u::Unitful.Units{<:Any,Unitful.𝐋}, pt::AbstractCartesianPoint{<:Quantity{<:Real, Unitful.𝐋}})
    CartesianPoint(uconvert(u, get_x(pt)), uconvert(u, get_y(pt)), uconvert(u, get_z(pt)))
end

function Unitful.uconvert(u::Unitful.Units{<:Any,Unitful.NoDims}, pt::AbstractCartesianPoint{<:Quantity{<:Real, Unitful.NoDims}})
    CartesianPoint(uconvert(u, get_x(pt)), uconvert(u, get_y(pt)), uconvert(u, get_z(pt)))
end

Unitful.ustrip(pt::AbstractCartesianPoint) = CartesianPoint(ustrip(get_x(pt)), ustrip(get_y(pt)), ustrip(get_z(pt)))






"""
    struct CartesianPoint{T} <: AbstractCartesianPoint{T}

Describes a three-dimensional point in Cartesian coordinates.

## Fields 
* `x`: x-coordinate (in m).
* `y`: y-coordinate (in m).
* `z`: z-coordinate (in m).

Given a point `pt = CartesianPoint(x, y, z)`, use `pt - cartesian_zero` to get
a `CartesianVector` from the origin to point `pt`.

See also [`CylindricalPoint`](@ref).
"""
struct CartesianPoint{T} <: AbstractCartesianPoint{T}
    x::T
    y::T
    z::T
end

#Type promotion happens here
function CartesianPoint(x::TX, y::TY, z::TZ) where {TX<:Real,TY<:Real,TZ<:Real}
    # ToDo: Simplify this:
    eltypes = _csg_get_promoted_eltype.((TX,TY,TZ))
    T = float(promote_type(eltypes...))
    CartesianPoint{T}(T(x),T(y),T(z))
end

CartesianPoint(; x = 0, y = 0, z = 0) = CartesianPoint(x, y, z)

CartesianPoint{T}(;x = 0, y = 0, z = 0) where {T} = CartesianPoint{T}(T(x),T(y),T(z))

# ToDo: Remove this, if possible
CartesianPoint{T}(v::AbstractVector) where {T} = cartesian_zero + v

CartesianPoint{T}(pt::CartesianPoint{T}) where {T} = pt

function CartesianPoint{T}(pt::CartesianPoint{U}) where {T, U} 
    return CartesianPoint{T}(convert(T, pt.x), convert(T, pt.y), convert(T, pt.z))
end

function Base.convert(::Type{CartesianPoint{T}}, pt::CartesianPoint) where {T}
    return CartesianPoint{T}(pt)
end


@inline get_x(pt::CartesianPoint) = pt.x
@inline get_y(pt::CartesianPoint) = pt.y
@inline get_z(pt::CartesianPoint) = pt.z


Base.keys(::CartesianPoint) = (:x, :y, :z)
Base.getindex(p::CartesianPoint, k::Symbol) = getfield(p, k)

AbstractCoordinatePoint{T, Cartesian}(x::Real, y::Real, z::Real) where T = CartesianPoint{T}(x, y, z)


Base.zero(::CartesianPoint{T}) where {T} = CartesianPoint{T}(zero(T),zero(T),zero(T))
# @inline Base.Tuple(pt::CartesianPoint) = (pt.x, pt.y, pt.z)
@inline Base.copy(pt::CartesianPoint) = CartesianPoint(pt.x, pt.y, pt.z)

@inline Base.:(+)(pt::CartesianPoint, v::CartesianVector) = CartesianPoint(pt.x + v.x, pt.y + v.y, pt.z + v.z)
@inline Base.:(-)(pt::CartesianPoint, v::CartesianVector) = CartesianPoint(pt.x - v.x, pt.y - v.y, pt.z - v.z)

@inline function Base.:(+)(pt::CartesianPoint, v::AbstractVector)
    length(v) == 3 || throw(DimensionMismatch("Can only add vectors of length 3 to Cartesian points."))
    CartesianPoint(pt.x + v[1], pt.y + v[2], pt.z + v[3])
end

@inline function Base.:(-)(pt::CartesianPoint, v::AbstractVector)
    length(v) == 3 || throw(DimensionMismatch("Can only add vectors of length 3 to Cartesian points."))
    CartesianPoint(pt.x - v[1], pt.y - v[2], pt.z - v[3])
end

@inline Base.:(-)(a::CartesianPoint, b::CartesianPoint) = CartesianVector(a.x - b.x, a.y - b.y, a.z - b.z)


# Barycentric combination
function  Statistics.mean(A::AbstractArray{<:CartesianPoint}; dims = :)
    cartesian_zero + mean(_vec_from_zero, A; dims = dims)
end

# Base.broadcastable(p::CartesianPoint) = (p.x, p.y, p.z)

_vec_from_zero(pt::CartesianPoint) = pt - cartesian_zero



"""
    CartesianZero{T} <: AbstractCoordinatePoint{T, Cartesian}

Represents origin of the Cartesian coordinate system.

See also [`cartesian_zero`](@ref) and [`CartesianVector`](@ref).

Constructors:

```julia
CartesianZero{T}()
CartesianZero() == CartesianZero{Bool}()
```
"""
struct CartesianZero{T} <: AbstractCoordinatePoint{T, Cartesian} end
CartesianZero() = CartesianZero{Bool}()

@inline get_x(pt::CartesianZero{T}) where T = zero(T)
@inline get_y(pt::CartesianZero{T}) where T = zero(T)
@inline get_z(pt::CartesianZero{T}) where T = zero(T)


"""
    const cartesian_zero = CartesianZero()

Origin of the Cartesian coordinate system.

See also [`CartesianZero`](@ref) and [`CartesianVector`](@ref).
"""
const cartesian_zero = CartesianZero()


@inline Base.:(==)(::CartesianZero, ::CartesianZero) = true
@inline Base.isapprox(::CartesianZero, ::CartesianZero; kwargs...) = true

@inline Base.:(-)(pt::CartesianPoint, ::CartesianZero) = CartesianVector(pt.x, pt.y, pt.z)
@inline Base.:(-)(::CartesianZero, pt::CartesianPoint) = CartesianVector(-pt.x, -pt.y, -pt.z)


@inline Base.:(+)(::CartesianZero, v::CartesianVector) = CartesianPoint(v.x, v.y, v.z)
@inline Base.:(-)(::CartesianZero, v::CartesianVector) = CartesianPoint(-v.x, -v.y, -v.z)

@inline Base.:(+)(::CartesianZero, v::StaticVector{3}) = CartesianPoint(v[1], v[2], v[3])
@inline Base.:(-)(::CartesianZero, v::StaticVector{3}) = CartesianPoint(-v[1], -v[2], -v[3])

@inline function Base.:(+)(::CartesianZero, v::AbstractVector)
    length(v) == 3 || throw(DimensionMismatch("Can only subtract vectors of length 3 from Cartesian points."))
    CartesianPoint(v[1], v[2], v[3])
end

@inline function Base.:(-)(::CartesianZero, v::AbstractVector)
    length(v) == 3 || throw(DimensionMismatch("Can only add vectors of length 3 to Cartesian points."))
    CartesianPoint(-v[1], -v[2], -v[3])
end

@inline function Base.:(-)(::CartesianZero{T}, ::CartesianZero{U}) where {T,U}
    R = promote_type(T, U)
    return CartesianVector(zero(R), zero(R), zero(R))
end

@inline Base.zero(::CartesianZero{T}) where {T} = CartesianZero{T}()
@inline Base.iszero(::CartesianZero) = true

Base.:(*)(::CartesianZero{T}, u::Unitful.Units{<:Any,Unitful.𝐋}) where {T<:Real} = CartesianZero{Quantity{T, Unitful.𝐋, typeof(u)}}()

function Unitful.uconvert(u::Unitful.Units{T,Unitful.NoDims}, pt::CartesianZero{<:Quantity{U, Unitful.NoDims}}) where {T,U}
    CartesianZero{promote_type(T,U)}()
end



"""
    const AbstractCylindricalPoint{T} = AbstractCoordinatePoint{T,Cylindrical}

Supertype for cylindrical point types.
"""
const AbstractCylindricalPoint{T} = AbstractCoordinatePoint{T,Cylindrical}



"""
    struct CylindricalPoint{T} <: AbstractCylindricalPoint{T}

Describes a three-dimensional point in cylindrical coordinates. 

## Fields
* `r`: Radius (in m).
* `φ`: Polar angle (in rad).
* `z`: `z`-coordinate (in m).

!!! note 
    `φ == 0` corresponds to the `x`-axis in the Cartesian coordinate system.
    
See also [`CartesianPoint`](@ref).
"""
struct CylindricalPoint{T} <: AbstractCylindricalPoint{T}
    r::T
    φ::T
    z::T
    CylindricalPoint{T}(r::T, φ::T, z::T) where {T<:AbstractFloat} = new(r, mod(φ,T(2π)), z)
    CylindricalPoint{T}(r::Real, φ::Real, z::Real) where {T} = new(T(r), mod(T(φ),T(2π)), T(z))
end

function CylindricalPoint(r::TR, φ::TP, z::TZ) where {TR<:Real,TP<:Real,TZ<:Real}
    # ToDo: Simplify this:
    eltypes = _csg_get_promoted_eltype.((TR,TP,TZ))
    T = float(promote_type(eltypes...))
    CylindricalPoint{T}(T(r),T(φ),T(z))
end

CylindricalPoint(; r = 0, φ = 0, z = 0) = CylindricalPoint(r,φ,z)

CylindricalPoint{T}(; r = 0, φ = 0, z = 0) where {T} = CylindricalPoint{T}(T(r),T(φ),T(z))

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

AbstractCoordinatePoint{T, Cylindrical}(r::Real, φ::Real, z::Real) where T = CylindricalPoint{T}(r, φ, z)

Base.:(==)(a::CylindricalPoint, b::CylindricalPoint) = a.r == b.r && a.φ == b.φ && a.z == b.z

function Base.isapprox(a::CylindricalPoint, b::CylindricalPoint; kwargs...)
    return isapprox(a.r, b.r; kwargs...) && isapprox(a.φ, b.φ; kwargs...) && isapprox(a.z, b.z; kwargs...)
end

Base.zero(::CylindricalPoint{T}) where {T} = CylindricalPoint{T}(zero(T), zero(T), zero(T))
Base.iszero(pt::CylindricalPoint) = iszero(pt.r) && iszero(pt.z)

@inline Base.copy(pt::CylindricalPoint) = CylindricalPoint(pt.r, pt.φ, pt.z)

@inline Base.:(+)(pt::CylindricalPoint, v::CartesianVector) = CylindricalPoint(CartesianPoint(pt) + v)

@inline Base.:(-)(pt::CylindricalPoint, v::CartesianVector) = CylindricalPoint(CartesianPoint(pt) - v)

@inline Base.:(-)(a::CylindricalPoint, b::CylindricalPoint) = CartesianPoint(a) - CartesianPoint(b)

# Barycentric combination
function  Statistics.mean(A::AbstractArray{<:CylindricalPoint}; dims = :)
    CylindricalPoint(cartesian_zero + mean(_vec_from_zero, A; dims = dims))
end

_vec_from_zero(pt::CylindricalPoint) = CartesianPoint(pt) - cartesian_zero


# function _Δφ(φ1::T, φ2::T)::T where {T}
#     δφ = mod(φ2 - φ1, T(2π))
#     min(δφ, T(2π) - δφ)
# end
# 
# _φNear(φ::Real, φMin::T, φMax::T) where {T} = _Δφ(T(φ),φMin) ≤ _Δφ(T(φ),φMax) ? φMin : φMax
