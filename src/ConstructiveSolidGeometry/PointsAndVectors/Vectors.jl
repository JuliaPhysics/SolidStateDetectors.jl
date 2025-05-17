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
"""
struct CartesianVector{T} <: AbstractCoordinateVector{T, Cartesian}
    x::T
    y::T
    z::T
end

#Type promotion happens here
function CartesianVector(x::TX, y::TY, z::TZ) where {TX<:Real,TY<:Real,TZ<:Real}
    eltypes = _csg_get_promoted_eltype.((TX,TY,TZ))
    T = float(promote_type(eltypes...))
    CartesianVector{T}(T(x),T(y),T(z))
end

function CartesianVector(;
    x = 0,
    y = 0,
    z = 0
)
    CartesianVector(x,y,z)
end

function CartesianVector{T}(;
    x = 0,
    y = 0,
    z = 0
) where {T}
    CartesianVector{T}(T(x),T(y),T(z))
end

zero(VT::Type{<:AbstractCoordinateVector{T}}) where {T} = VT(zero(T),zero(T),zero(T))

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
