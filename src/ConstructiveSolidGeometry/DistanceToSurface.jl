function distance_to_surface(point::AbstractCoordinatePoint{T}, c::Union{<:AbstractConstructiveGeometry{T}, <:AbstractVolumePrimitive{T}}) where {T}
    min_distance::T = T(Inf)
    for surf in surfaces(c)
        d::T = distance_to_surface(point, surf)
        min_distance = ifelse(isnan(d), min_distance, min(d, min_distance))
    end
    @assert min_distance!=Inf "Cound determine the distance to surface from this point. Something unexpected issue happens with the geometry or this point."
    return min_distance
end