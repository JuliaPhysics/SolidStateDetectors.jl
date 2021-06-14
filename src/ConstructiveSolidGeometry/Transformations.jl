const Transformations{T} = NamedTuple{(:rotation, :translation), Tuple{SMatrix{3, 3, T, 9}, CartesianVector{T}}}
# Transformation rule:
#   1. Rotate
#   2. Translate
# Not sure yet if Rotations will always be of type SMatrix{3, 3, T, 9} -> Any

Transformations{T}() where {T} = (rotation = one(SMatrix{3, 3, T, 9}), translation = zero(CartesianVector{T}))

rotate(p::VP, r::AbstractMatrix) where {VP <: AbstractVolumePrimitive} = VP(p, origin = r * p.origin, rotation = r * p.rotation)
(*)(r::AbstractMatrix, p::AbstractVolumePrimitive) = rotate(p, r)

translate(p::VP, v::CartesianVector) where {VP <: AbstractVolumePrimitive} = VP(p, origin = p.origin + v, rotation = p.rotation)
(+)(p::AbstractVolumePrimitive, v::CartesianVector) = translate(p, v)

transform(g::AbstractPrimitive, t::Transformations) =
    translate(rotate(g, t.rotation), t.translation)

transform(csg::G, t::Transformations) where {G <: AbstractConstructiveGeometry} = G(transform(csg.a, t), transform(csg.b, t))


function combine_transformations(current::Transformations, new::Transformations)
    r = new.rotation * current.rotation    
    t = new.rotation * current.translation
    t = t + new.translation
    (rotation = r, translation = t)
end

