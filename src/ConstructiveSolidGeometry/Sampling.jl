function sample_surface(c::AbstractVolumePrimitive{T}, sampling) where {T}
    samples = [ point
    for surf in get_decomposed_surfaces(c)
    for point in sample(surf, sampling)  ]
end

function sample_surface(c::AbstractSurfacePrimitive{T}, sampling) where {T}
    sample(c, sampling)
end

sample(sg::ScaledGeometry{T}, sampling) where {T} = scale!(sample(sg.p, sampling), inv.(sg.inv_s))
sample(rg::RotatedGeometry{T}, sampling) where {T} = rotate!(sample(rg.p, sampling), inv(rg.inv_r))
sample(tg::TranslatedGeometry{T}, sampling) where {T} = translate!(sample(tg.p, sampling), tg.t)
sample(csg::AbstractConstructiveGeometry{T}, sampling) where {T} = vcat(sample(csg.a, sampling), sample(csg.b, sampling))

sample_surface(sg::ScaledGeometry{T}, sampling) where {T} = scale!(sample_surface(sg.p, sampling), inv.(sg.inv_s))
sample_surface(rg::RotatedGeometry{T}, sampling) where {T} = rotate!(sample_surface(rg.p, sampling), inv(rg.inv_r))
sample_surface(tg::TranslatedGeometry{T}, sampling) where {T} = translate!(sample_surface(tg.p, sampling), tg.t)
sample_surface(csg::AbstractConstructiveGeometry{T}, sampling) where {T} = vcat(sample_surface(csg.a, sampling), sample_surface(csg.b, sampling))