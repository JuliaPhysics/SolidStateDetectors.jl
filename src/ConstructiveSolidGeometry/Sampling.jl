function sample_surface(c::AbstractVolumePrimitive{T}, sampling) where {T}
    samples = [ point
    for surf in get_decomposed_surfaces(c)
    for point in sample(surf, sampling)  ]
end

function sample_surface(c::AbstractSurfacePrimitive{T}, sampling) where {T}
    sample(c, sampling)
end
