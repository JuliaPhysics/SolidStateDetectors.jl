# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

abstract type AbstractCoordinateVector{T, N, S} <: StaticArrays.FieldVector{N, T} end

struct CartesianVector{ T <: RealQuantity } <: AbstractCoordinateVector{T, 3, :cartesian}
    x::T
    y::T
    z::T
end

function (+)(p::CartesianPoint{T}, v::CartesianVector{T})::CartesianPoint{T} where {T <: Real}
    return CartesianPoint{T}( p.x + v.x, p.y + v.y, p.z + v.z )
end
function (+)(u::CartesianVector{T}, v::CartesianVector{T})::CartesianVector{T} where {T <: Real}
    return CartesianVector{T}( u.x + v.x, u.y + v.y, u.z + v.z )
end
function broadcast(f::Function, v::Vector{CartesianPoint{T}}, p::CartesianVector{T}) where {T}
    f(v ,[p for i in eachindex(v)])
end

function (-)(p::CartesianPoint{T}, v::CartesianVector{T})::CartesianPoint{T} where {T <: Real}
    return CartesianPoint{T}( p.x - v.x, p.y - v.y, p.z - v.z )
end
function (-)(u::CartesianVector{T}, v::CartesianVector{T})::CartesianVector{T} where {T <: Real}
    return CartesianVector{T}( u.x - v.x, u.y - v.y, u.z - v.z )
end
function (-)(p1::CartesianPoint{T}, p2::CartesianPoint{T})::CartesianVector{T} where {T <: Real}
    return CartesianVector{T}( p1.x - p2.x, p1.y - p2.y, p1.z - p2.z )
end


struct CylindricalVector{ T <: RealQuantity } <: AbstractCoordinateVector{T, 3, :cylindrical}
    r::T
    φ::T # in radian
    z::T
end

function (+)(p::CylindricalPoint{T}, v::CylindricalVector{T})::CylindricalPoint{T} where {T <: Real}
    error("Not yet defined")
end
function (-)(p::CylindricalPoint{T}, v::CylindricalVector{T})::CylindricalPoint{T} where {T <: Real}
    error("Not yet defined")
end
function (-)(p1::CylindricalPoint{T}, p2::CylindricalPoint{T})::CylindricalPoint{T} where {T <: Real}
    error("Not yet defined")
end