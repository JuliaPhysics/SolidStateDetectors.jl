# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

abstract type AbstractCoordinatePoint{T, N, S} <: StaticArrays.FieldVector{N, T} end

struct CartesianPoint{ T <: RealQuantity } <: AbstractCoordinatePoint{T, 3, :cartesian}
    x::T
    y::T
    z::T
end

struct CylindricalPoint{ T <: RealQuantity } <: AbstractCoordinatePoint{T, 3, :cylindrical}
    r::T
    φ::T # in radian
    z::T
    CylindricalPoint{T}(r,φ,z) where {T <: RealQuantity} = new(r,mod(φ,2π),z)
end


function CylindricalPoint( pt::CartesianPoint{T} )::CylindricalPoint{T} where {T <: RealQuantity}
    return CylindricalPoint{T}(sqrt(pt.x * pt.x + pt.y * pt.y), atan(pt.y, pt.x), pt.z)
end

function CartesianPoint( pt::CylindricalPoint{T} )::CartesianPoint{T} where {T <: RealQuantity}
    sφ::T, cφ::T = sincos(pt.φ)
    return CartesianPoint{T}(pt.r * cφ, pt.r * sφ, pt.z)
end

function convert(type::Type{CylindricalPoint}, origin::SolidStateDetectors.CartesianPoint{T})::CylindricalPoint{T} where T
    return CylindricalPoint(origin)
end

function convert(type::Type{CartesianPoint}, origin::SolidStateDetectors.CylindricalPoint{T})::CartesianPoint{T} where T
    return CartesianPoint(origin)
end

function (+)(p1::CartesianPoint{T}, p2::CartesianPoint{T})::CartesianPoint{T} where {T <: Real}
    return CartesianPoint{T}( p1.x + p2.x, p1.y + p2.y, p1.z + p2.z )
end

# function (-)(p1::CartesianPoint{T}, p2::CartesianPoint{T})::CartesianPoint{T} where {T <: Real}
#     return CartesianPoint{T}( p1.x - p2.x, p1.y - p2.y, p1.z - p2.z )
# end

function (+)(p1::CylindricalPoint{T}, p2::CylindricalPoint{T})::CylindricalPoint{T} where {T <: Real}
    return CylindricalPoint(CartesianPoint(p1) + CartesianPoint(p2))
end

function broadcast(f::Function, v::Vector{CartesianPoint{T}}, p::CartesianPoint{T}) where {T}
    f(v ,[p for i in eachindex(v)])
end

function broadcast(f::Function, v::Vector{CylindricalPoint{T}}, p::CylindricalPoint{T}) where {T}
    f(v ,[p for i in eachindex(v)])
end

@recipe function f(p::CylindricalPoint)
    @series begin
        CartesianPoint(p)
    end
end
@recipe function f(p::CartesianPoint)
    # st -> :scatter
    @series begin
        [p.x], [p.y], [p.z]
    end
end
@recipe function f(v::AbstractVector{<:CartesianPoint})
    st -> :scatter
    @series begin
        [v[i].x for i in eachindex(v) ], [v[i].y for i in eachindex(v) ], [v[i].z for i in eachindex(v) ]
    end
end
@recipe function f(v::AbstractVector{<:CylindricalPoint})
    st -> :scatter
    @series begin
        map(x->CartesianPoint(x),v)
    end
end
 ##aliases for ease of use
const cyp = CylindricalPoint
const cap = CartesianPoint
