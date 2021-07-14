sample(csg::AbstractConstructiveGeometry{T}, sampling...) where {T} = vcat(sample(csg.a, sampling...), sample(csg.b, sampling...))

# for sampling surfaces with predefined ticks
get_r_ticks(a::AbstractSurfacePrimitive{T}, g::CylindricalTicksTuple{T}) where {T} = _get_ticks(g.r, get_r_limits(a)...)
get_φ_ticks(a::AbstractSurfacePrimitive{T}, g::CylindricalTicksTuple{T}) where {T} = _get_ticks(g.φ, get_φ_limits(a)...)
get_z_ticks(a::AbstractSurfacePrimitive{T}, g::CylindricalTicksTuple{T}) where {T} = _get_ticks(g.z, get_z_limits(a)...)
get_z_ticks(a::AbstractSurfacePrimitive{T}, g::CartesianTicksTuple{T}) where {T} = _get_ticks(g.z, get_z_limits(a)...)
_get_ticks(ticks::Vector{T}, Min::T, Max::T, args...) where {T} = (ticks[1] <= Max && ticks[end] >= Min) ? vcat(Min, filter(tick -> Min < tick < Max, ticks), Max) : []
