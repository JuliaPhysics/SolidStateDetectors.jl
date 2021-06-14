const Transformations{T} = NamedTuple{(:scale, :rotation, :translation), 
    Tuple{SVector{3, T}, SMatrix{3, 3, T, 9}, CartesianVector{T}}}
# Transformation rule:
#   1. Scale
#   2. Rotate
#   3. Translate

scale(p::VP, s::SVector{3, <:Any}) where {VP <: AbstractVolumePrimitive} = VP(p, scaling = s, origin = p.origin, rotation = p.rotation)
(*)(s::SVector{3, <:Any}, p::AbstractVolumePrimitive) = scale(p, s)

rotate(p::VP, r::AbstractMatrix) where {VP <: AbstractVolumePrimitive} = VP(p, origin = r * p.origin, rotation = r * p.rotation)
(*)(r::AbstractMatrix, p::AbstractVolumePrimitive) = rotate(p, r)

translate(p::VP, v::CartesianVector) where {VP <: AbstractVolumePrimitive} = VP(p, origin = p.origin + v, rotation = p.rotation)
(+)(p::AbstractVolumePrimitive, v::CartesianVector) = translate(p, v)

transform(g::AbstractPrimitive, t::Transformations) =
    translate(rotate(scale(g, t.scale), t.rotation), t.translation)

transform(csg::G, t::Transformations) where {G <: AbstractConstructiveGeometry} = G(transform(csg.a, t), transform(csg.b, t))

# Not sure yet if Rotations will always be of type SMatrix{3, 3, T, 9} -> Any

function combine_transformations(current::Transformations, new::Transformations)
    s = current.scale .* new.scale
    t = scale(current.translation, new.scale)
    r = new.rotation * current.rotation    
    t = new.rotation * t
    t = t + new.translation
    (scale = s, rotation = r, translation = t)
end

