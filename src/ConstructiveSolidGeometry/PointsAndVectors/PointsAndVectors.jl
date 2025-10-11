abstract type AbstractCoordinatePoint{T, S} <: StaticArrays.FieldVector{3, T} end
abstract type AbstractCoordinateVector{T, S} <: StaticArrays.FieldVector{3, T} end

include("Points.jl")
include("Vectors.jl")

@inline distance_squared(v::CartesianVector) = v.x * v.x + v.y * v.y + v.z * v.z
@inline distance_squared(p1::CartesianPoint{T}, p2::CartesianPoint{T}) where {T <: Real} = distance_squared(CartesianVector(p1 .- p2))
export distance_squared