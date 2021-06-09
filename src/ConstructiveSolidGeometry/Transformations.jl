transform(csg::G, t) where {G <: AbstractConstructiveGeometry} = G(transform(csg.a, t), transform(csg.b, t))
transform(csg::G, t::Missing) where {G <: AbstractConstructiveGeometry} = csg

transform(g::AbstractPrimitive, v::AbstractVector) =
    length(v) > 1 ? transform(transform(g, v[end]), v[1:end-1]) : transform(g, v[1])

transform(p::AbstractVolumePrimitive, v::CartesianVector) = p + v
transform(p::AbstractVolumePrimitive, r::AbstractMatrix) = r * p
transform(p::AbstractVolumePrimitive, s::SVector{3}) = s * p
transform(p::AbstractVolumePrimitive, ::Missing) = p