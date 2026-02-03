function distance_to_surface(point::AbstractCoordinatePoint{T}, c::Union{<:AbstractConstructiveGeometry{T}, <:AbstractVolumePrimitive{T}}) where {T}
    # determine the minimum distance to all surfaces of the geometry
    distances = broadcast(surf -> let d = distance_to_surface(point, surf); isfinite(d) ? d : T(Inf) end, surfaces(c))
    min_distance::T = min(distances..., T(Inf))
    isfinite(min_distance) || throw(AssertionError("Couldn't determine the distance to surface from this point. Some unexpected issue happens with the geometry or this point."))
    return min_distance
end
