sample_surface(c::AbstractVolumePrimitive{T}, sampling...) where {T} = [ point for surf in get_decomposed_surfaces(c) for point in sample(surf, sampling...)  ]
sample_surface(c::AbstractSurfacePrimitive{T}, sampling...) where {T} = sample(c, sampling...)

sample(sg::ScaledGeometry{T}, sampling...) where {T} = scale!(sample(sg.p, sampling...), inv.(sg.inv_s))
sample(rg::RotatedGeometry{T}, sampling...) where {T} = rotate!(sample(rg.p, sampling...), inv(rg.inv_r))
sample(tg::TranslatedGeometry{T}, sampling...) where {T} = translate!(sample(tg.p, sampling...), tg.t)
sample(csg::AbstractConstructiveGeometry{T}, sampling...) where {T} = vcat(sample(csg.a, sampling...), sample(csg.b, sampling...))
sample(tg::TranslatedGeometry{T}, g::CartesianTuple{T}) where {T} = translate!(sample(tg.p, (x = g.x .- tg.t[1], y = g.y .- tg.t[2], z = g.z .- tg.t[3])), tg.t)
sample(tg::TranslatedGeometry{T}, g::CylindricalTuple{T}) where {T} = sample(tg.p, (r = g.r, φ = g.φ, z = g.z .- tg.t[3])) .+ (CartesianVector{T}(0,0,tg.t[3]),)

sample_surface(sg::ScaledGeometry{T}, sampling...) where {T} = scale!(sample_surface(sg.p, sampling...), inv.(sg.inv_s))
sample_surface(rg::RotatedGeometry{T}, sampling...) where {T} = rotate!(sample_surface(rg.p, sampling...), inv(rg.inv_r))
sample_surface(tg::TranslatedGeometry{T}, sampling...) where {T} = translate!(sample_surface(tg.p, sampling...), tg.t)
sample_surface(csg::AbstractConstructiveGeometry{T}, sampling...) where {T} = vcat(sample_surface(csg.a, sampling...), sample_surface(csg.b, sampling...))
sample_surface(tg::TranslatedGeometry{T}, g::CartesianTuple{T}) where {T} = translate!(sample_surface(tg.p, (g[1] .- tg.t[1], g[2] .- tg.t[2], g[3] .- tg.t[3])), tg.t)
sample_surface(tg::TranslatedGeometry{T}, g::CylindricalTuple{T}) where {T} = sample_surface(tg.p, (g[1], g[2], g[3] .- tg.t[3])) .+ (CartesianVector{T}(0,0,tg.t[3]),)

# for sampling surfaces with predefined ticks
get_r_ticks(a::AbstractSurfacePrimitive{T}, g::CylindricalTuple{T}) where {T} = _get_ticks(g.r, get_r_limits(a)...)
get_φ_ticks(a::AbstractSurfacePrimitive{T}, g::CylindricalTuple{T}) where {T} = _get_ticks(g.φ, get_φ_limits(a)...)
get_z_ticks(a::AbstractSurfacePrimitive{T}, g::CylindricalTuple{T}) where {T} = _get_ticks(g.z, get_z_limits(a)...)
get_z_ticks(a::AbstractSurfacePrimitive{T}, g::CartesianTuple{T}) where {T} = _get_ticks(g.z, get_z_limits(a)...)
_get_ticks(ticks::Vector{T}, Min::T, Max::T, args...) where {T} = (ticks[1] <= Max && ticks[end] >= Min) ? vcat(Min, filter(tick -> Min < tick < Max, ticks), Max) : []
