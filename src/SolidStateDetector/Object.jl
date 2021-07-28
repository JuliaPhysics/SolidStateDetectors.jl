abstract type AbstractObject{T <: SSDFloat} end

@inline function in(pt::AbstractCoordinatePoint{T}, c::AbstractObject{T})::Bool where {T}
    return in(pt, c.geometry)
end

function in(pt::AbstractCoordinatePoint{T}, v::AbstractVector{<:AbstractObject{T}}) where {T <: SSDFloat}
    reduce((x, object) -> x || in(pt, object), v, init = false)
end
