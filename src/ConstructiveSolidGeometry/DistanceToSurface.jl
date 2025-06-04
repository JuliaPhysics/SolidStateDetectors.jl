function distance_to_surface(point::AbstractCoordinatePoint{T}, c::AbstractGeometry{T}) where {T}
    distance_list=[distance_to_surface(point,surf) for surf in surfaces(c)]
    minimum(filter(!isnan,distance_list))
end
