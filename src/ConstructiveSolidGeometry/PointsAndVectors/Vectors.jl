# struct PlanarVector{T} <: AbstractPlanarVector{T}
#     u::T
#     v::T
# end

"""
    struct CartesianVector{T} <: AbstractCoordinateVector{T, Cartesian}

Describes a three-dimensional vector in Cartesian coordinates.

## Fields 
* `x`: x-coordinate (in m).
* `y`: y-coordinate (in m).
* `z`: z-coordinate (in m).

Given a vector `v = CartesianPoint(Δx, Δy, Δz)`, use `cartesian_zero + v` to
get the `CartesianPoint(0 + Δx, 0 + Δy, 0 + Δz)`.

See also [`CartesianPoint`](@ref) .
"""
struct CartesianVector{T} <: AbstractCoordinateVector{T, Cartesian}
    x::T
    y::T
    z::T

    CartesianVector{T}(x::T, y::T, z::T) where {T<:AbstractFloat} = new(x, y, z)
    CartesianVector{T}(x::Real, y::Real, z::Real) where {T} = new(T(x), T(y), T(z))
end

function CartesianVector(x, y, z)
    for (name, val) in zip((:x, :y, :z), (x, y, z))
        (isa(val, Real) || isa(val, Unitful.Length)) ||
            throw(ArgumentError("Expected $(name) to be a Real or Unitful.Quantity, got $(typeof(val))"))
    end
    x_val = to_internal_units(x)
    y_val = to_internal_units(y)
    z_val = to_internal_units(z)

    return CartesianVector(x_val, y_val, z_val)
end

#Type promotion happens here
function CartesianVector(x::TX, y::TY, z::TZ) where {TX<:Real,TY<:Real,TZ<:Real}
    # ToDo: Simplify this:
    eltypes = _csg_get_promoted_eltype.((TX,TY,TZ))
    T = float(promote_type(eltypes...))
    CartesianVector{T}(T(x),T(y),T(z))
end

CartesianVector(; x = 0, y = 0, z = 0) = CartesianVector(x,y,z)

CartesianVector{T}(; x = 0, y = 0, z = 0) where {T} = CartesianVector{T}(T(x),T(y),T(z))

@inline Base.zero(::Union{<:CartesianVector{T}, <:Type{CartesianVector{T}}}) where {T} = CartesianVector{T}(zero(T),zero(T),zero(T))
@inline Base.iszero(v::CartesianVector) = iszero(v.x) && iszero(v.y) && iszero(v.z)

# Need to specialize multiplication and division with units for CartesianVector, otherwise result would just be an SVector:
Base.:(*)(v::CartesianVector{<:Real}, u::Unitful.Units{<:Any,Unitful.𝐋}) = CartesianVector(v.x * u, v.y * u, v.z * u)
#Base.:(/)(v::CartesianVector{<:Quantity{<:Real, Unitful.𝐋}}, u::Unitful.Units) = CartesianVector(v.x / u, v.y / u, v.z / u)

# For user convenience, to enable constructs like `ustrip(u"mm", v)` and `NoUnits(v / u"mm")`, we'll support uconvert/ustrip
# for CartesianVector. Unitful doesn't encourage defining those for collections, but we'll view cartesian vectors as single
# mathematical objects in regard to units:

#ToDo: Uncomment when units are supported
#Unitful.uconvert(u::Unitful.Units{<:Any,Unitful.𝐋}, v::CartesianVector{<:Quantity{<:Real, Unitful.𝐋}}) = CartesianVector(uconvert(u, v.x), uconvert(u, v.y), uconvert(u, v.z))
#Unitful.uconvert(u::Unitful.Units{<:Any,Unitful.NoDims}, v::CartesianVector{<:Quantity{<:Real, Unitful.NoDims}}) = CartesianVector(uconvert(u, v.x), uconvert(u, v.y), uconvert(u, v.z))
#Unitful.ustrip(v::CartesianVector) = CartesianVector(ustrip(v.x), ustrip(v.y), ustrip(v.z))


# @inline rotate(pt::CartesianPoint{T}, r::RotMatrix{3,T,TT}) where {T, TT} = r.mat * pt
# @inline rotate(pt::CylindricalPoint{T}, r::RotMatrix{3,T,TT}) where {T, TT} = CylindricalPoint(rotate(CartesianPoint(pt), r))
# @inline rotate!(vpt::Vector{<:AbstractCoordinatePoint{T}}, r::RotMatrix{3,T,TT}) where {T, TT} = begin for i in eachindex(vpt) vpt[i] = rotate(vpt[i], r) end; vpt end
# @inline rotate!(vvpt::Vector{<:Vector{<:AbstractCoordinatePoint{T}}}, r::RotMatrix{3,T,TT}) where {T, TT} = begin for i in eachindex(vvpt) rotate!(vvpt[i], r) end; vvpt end

# @inline scale(pt::CartesianPoint{T}, v::SVector{3,T}) where {T} = pt .* v
# @inline scale(pt::CylindricalPoint{T}, v::SVector{3,T}) where {T} = CylindricalPoint(scale(CartesianPoint(pt),v))
# @inline scale!(vpt::Vector{<:AbstractCoordinatePoint{T}}, v::SVector{3,T}) where {T} = begin for i in eachindex(vpt) vpt[i] = scale(vpt[i], v) end; vpt end
# @inline scale!(vvpt::Vector{<:Vector{<:AbstractCoordinatePoint{T}}}, v::SVector{3,T}) where {T} = begin for i in eachindex(vvpt) scale!(vvpt[i], v) end; vvpt end

# @inline translate(pt::CartesianPoint{T}, v::CartesianVector{T}) where {T} = pt + v
# @inline translate(pt::CylindricalPoint{T}, v::CartesianVector{T}) where {T} = CylindricalPoint(translate(CartesianPoint(pt), v))
# @inline translate!(vpt::Vector{<:AbstractCoordinatePoint{T}}, v::CartesianVector{T}) where {T} = begin for i in eachindex(vpt) vpt[i] = translate(vpt[i], v) end; vpt end
# @inline translate!(vvpt::Vector{<:Vector{<:AbstractCoordinatePoint{T}}}, v::CartesianVector{T}) where {T} =  begin for i in eachindex(vvpt) translate!(vvpt[i], v) end; vvpt end
