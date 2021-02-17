function distance_to_surface(point::AbstractCoordinatePoint{T}, c::AbstractVolumePrimitive{T}) where {T}
    minimum([ distance_to_surface(point,surf) for surf in get_decomposed_surfaces(c)])
end
