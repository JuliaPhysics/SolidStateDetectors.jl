abstract type AbstractCoordinatePoint{T, S} end

@inline Base.zero(PT::Type{<:AbstractCoordinatePoint{T}}) where {T} = PT(zero(T),zero(T),zero(T))

@inline StaticArrays.Size(::PT) where {PT<:AbstractCoordinatePoint} = Size(PT)
@inline StaticArrays.Size(::Type{PT}) where {PT<:AbstractCoordinatePoint} = Size{(fieldcount(PT),)}()
@inline Base.length(pt::AbstractCoordinatePoint) = prod(Size(pt))::Int
@inline Base.size(pt::AbstractCoordinatePoint) = Tuple(Size(pt))

Base.@propagate_inbounds Base.getindex(pt::AbstractCoordinatePoint, i::Int) = getfield(pt, i)


abstract type AbstractCoordinateVector{T, S} <: StaticArrays.FieldVector{3, T} end

include("Vectors.jl")
include("Points.jl")
include("AffineFrames.jl")

@inline distance_squared(v::CartesianVector) = v.x * v.x + v.y * v.y + v.z * v.z
@inline distance_squared(p1::CartesianPoint{T}, p2::CartesianPoint{T}) where {T <: Real} = distance_squared(p1 - p2)
