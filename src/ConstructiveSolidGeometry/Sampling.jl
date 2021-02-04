function sample_surface(c::AbstractVolumePrimitive{T}, sampling::Union{Real, NTuple{3,Int}}) where {T}
    samples = CylindricalPoint{T}[]
    for surf in get_decomposed_surfaces(c)
        append!(samples, sample(surf, sampling))
    end
    samples
end

function sample_surface(c::AbstractSurfacePrimitive{T}, sampling::Union{Real, NTuple{3,Int}}) where {T}
    sample(c, sampling)
end
