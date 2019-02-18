# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

abstract type CoordinateVector{T, N, S} <: StaticArrays.FieldVector{N, T} end

struct CartesianVector{ T <: RealQuantity } <: CoordinateVector{T, 3, :Cartesian}
    x::T
    y::T
    z::T
end

function (+)(p::CartesianPoint{T}, v::CartesianVector{T})::CartesianPoint{T} where {T <: Real} 
    return CartesianPoint{T}( p.x + v.x, p.y + v.y, p.z + v.z )
end
function (-)(p::CartesianPoint{T}, v::CartesianVector{T})::CartesianPoint{T} where {T <: Real} 
    return CartesianPoint{T}( p.x - v.x, p.y - v.y, p.z - v.z )
end
function (-)(p1::CartesianPoint{T}, p2::CartesianPoint{T})::CartesianVector{T} where {T <: Real} 
    return CartesianVector{T}( p1.x - p2.x, p1.y - p2.y, p1.z - p2.z )
end


struct CylindricalVector{ T <: RealQuantity } <: CoordinateVector{T, 3, :Cylindrical}
    r::T
    φ::T # in radian
    z::T
end


