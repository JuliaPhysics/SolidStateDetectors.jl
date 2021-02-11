function sample_surface(c::AbstractVolumePrimitive{T}, sampling) where {T}
    samples = AbstractCoordinatePoint{T}[]
    for surf in get_decomposed_surfaces(c)
        append!(samples, sample(surf, sampling))
    end
    samples
end

function sample_surface(c::AbstractSurfacePrimitive{T}, sampling) where {T}
    sample(c, sampling)
end
