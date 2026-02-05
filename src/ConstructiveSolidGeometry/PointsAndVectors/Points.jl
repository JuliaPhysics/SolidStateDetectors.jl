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


# Unitful uses the `*` and `/` operators to combine values with units. But mathematically that's really not an algebraic
# product (not defined for affine points), but a cartesian product, so supporting this should be fine:
Base.:(*)(pt::AbstractCartesianPoint{<:Real}, u::Unitful.Units{<:Any,Unitful.𝐋}) = CartesianPoint(get_x(pt) * u, get_y(pt) * u, get_z(pt) * u)
# ToDo: Uncomment once units in points are supported
#=
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
=#


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
    CartesianPoint{T}(x::T, y::T, z::T) where {T<:AbstractFloat} = new(x, y, z)
    CartesianPoint{T}(x::Real, y::Real, z::Real) where {T} = new(T(x), T(y), T(z))
end

#Units support
function CartesianPoint(x, y, z)
    for (name, pt) in zip((:x, :y, :z), (x, y, z))
        (pt isa Real || pt isa Unitful.Length) ||
            throw(ArgumentError("Expected $(name) to be a length or Real, got unit $(Unitful.unit(pt))"))
    end
    x_val = to_internal_units(x)
    y_val = to_internal_units(y)
    z_val = to_internal_units(z)

    return CartesianPoint(x_val, y_val, z_val)
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


Base.eltype(::Type{<:CartesianPoint{T}}) where {T} = T
Base.keys(::CartesianPoint) = (:x, :y, :z)
Base.getindex(p::CartesianPoint, k::Symbol) = getfield(p, k)

AbstractCoordinatePoint{T, Cartesian}(x::Real, y::Real, z::Real) where T = CartesianPoint{T}(x, y, z)


@inline Base.zero(::CartesianPoint{T}) where {T} = CartesianPoint{T}(zero(T),zero(T),zero(T))
@inline Base.iszero(pt::CartesianPoint) = iszero(pt.x) && iszero(pt.y) && iszero(pt.z)
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

@inline Base.transpose(pt::CartesianPoint) = CartesianPoint(transpose(pt.x), transpose(pt.y), transpose(pt.z))
@inline Base.adjoint(pt::CartesianPoint) = CartesianPoint(adjoint(pt.x), adjoint(pt.y), adjoint(pt.z))

barycenter(X, args...; kwargs...) = _barycenter_impl(eltype(X), X, args...; kwargs...)
function _barycenter_impl(::Type{T}, points; kwargs...) where {T<:AbstractCoordinatePoint}
    mean_vector = mean(Base.Fix2(-, cartesian_zero), CartesianPoint.(points); kwargs...)
    return convert(T, cartesian_zero + mean_vector)::T
end

function _barycenter_impl(::Type{T}, points, weights::StatsBase.AbstractWeights; kwargs...) where {T<:AbstractCoordinatePoint}
    mean_vector = mean(CartesianPoint.(points) .- Ref(cartesian_zero), weights; kwargs...)
    return convert(T, cartesian_zero + mean_vector)::T
end

# Fallback if a Vector of Vectors (e.g. StaticVectors) is passed to barycenter
function _barycenter_impl(::Type{T}, args...; kwargs...) where {T <: AbstractVector}
    convert(T, mean(args...; kwargs...))
end


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

@inline Base.zero(::Union{<:CartesianZero{T}, <:Type{CartesianZero{T}}}) where {T} = CartesianZero{T}()
@inline Base.iszero(::CartesianZero) = true


# ToDo: Revert this once we support units internally.
Base.:(*)(::CartesianZero{T}, u::Unitful.Units{<:Any,Unitful.𝐋}) where {T<:Real} = CartesianZero{Quantity{T, Unitful.𝐋, typeof(u)}}()
Base.:(*)(z::CartesianZero, u::Unitful.Quantity) = z

# ToDo: Uncomment once units in points are supported
#=
function Unitful.uconvert(u::Unitful.Units{T,Unitful.NoDims}, pt::CartesianZero{<:Quantity{U, Unitful.NoDims}}) where {T,U}
    CartesianZero{promote_type(T,U)}()
end
=#


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

#Units support
function CylindricalPoint(r, φ, z)
    if !(r isa Real || r isa Unitful.Length)
        throw(ArgumentError("Expected `r` to be a length or Real, got unit $(Unitful.unit(r))"))
    end

    if !(φ isa Real || φ isa Unitful.Quantity{<:Real, NoDims})
        throw(ArgumentError("Expected `φ` to be an angle or Real, got unit $(Unitful.unit(φ))"))
    end

    if !(z isa Real || z isa Unitful.Length)
        throw(ArgumentError("Expected `z` to be a length or Real, got unit $(Unitful.unit(z))"))
    end

    r_val = to_internal_units(r)
    φ_val = to_internal_units(φ)
    z_val = to_internal_units(z)

    return CylindricalPoint(r_val, φ_val, z_val)
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

function Base.convert(::Type{CylindricalPoint{T}}, pt::AbstractCoordinatePoint)::CylindricalPoint{T} where {T}
    return CylindricalPoint(pt)
end

function Base.convert(::Type{CartesianPoint{T}}, pt::AbstractCoordinatePoint)::CartesianPoint{T} where {T}
    return CartesianPoint(pt)
end

AbstractCoordinatePoint{T, Cylindrical}(r::Real, φ::Real, z::Real) where T = CylindricalPoint{T}(r, φ, z)

Base.:(==)(a::CylindricalPoint, b::CylindricalPoint) = a.r == b.r && a.φ == b.φ && a.z == b.z

function Base.isapprox(a::CylindricalPoint, b::CylindricalPoint; kwargs...)
    return isapprox(a.r, b.r; kwargs...) && isapprox(a.φ, b.φ; kwargs...) && isapprox(a.z, b.z; kwargs...)
end

@inline Base.zero(::CylindricalPoint{T}) where {T} = CylindricalPoint{T}(zero(T), zero(T), zero(T))
@inline Base.iszero(pt::CylindricalPoint) = iszero(pt.r) && iszero(pt.z)

@inline Base.copy(pt::CylindricalPoint) = CylindricalPoint(pt.r, pt.φ, pt.z)

@inline Base.:(+)(pt::CylindricalPoint, v::CartesianVector) = CylindricalPoint(CartesianPoint(pt) + v)

@inline Base.:(-)(pt::CylindricalPoint, v::CartesianVector) = CylindricalPoint(CartesianPoint(pt) - v)

@inline Base.:(-)(a::CylindricalPoint, b::CylindricalPoint) = CartesianPoint(a) - CartesianPoint(b)

@inline Base.transpose(pt::CylindricalPoint) = CylindricalPoint(transpose(pt.r), transpose(pt.φ), transpose(pt.z)) 

@inline Base.adjoint(pt::CylindricalPoint) = CylindricalPoint(adjoint(pt.r), adjoint(pt.φ), adjoint(pt.z))

function to_internal_units(pt::CartesianPoint)
    CartesianPoint(to_internal_units(pt.x), to_internal_units(pt.y), to_internal_units(pt.z))
end

function to_internal_units(pt::CylindricalPoint)
    _r = to_internal_units(pt.r)
    _φ = to_internal_units(pt.φ)
    _z = to_internal_units(pt.z)

    return CartesianPoint(_r * cos(_φ), _r * sin(_φ), _z)
end

function to_internal_units(pt::AbstractCoordinatePoint)
    throw(ArgumentError("Unsupported point type $(typeof(pt)). Expected CartesianPoint or CylindricalPoint."))
end
