abstract type AbstractObject{T <: SSDFloat} end

@inline function in(pt::AbstractCoordinatePoint{T}, c::AbstractObject{T})::Bool where {T}
    return in(pt, c.geometry)
end

function in(p::AbstractCoordinatePoint{T}, v::AbstractVector{<:AbstractObject{T}}) where {T <: SSDFloat}
    reduce((x, object) -> x || in(p, object), v, init = false)
end
