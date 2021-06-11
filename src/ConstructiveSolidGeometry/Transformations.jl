const Transformations{T} = NamedTuple{(:scale, :rotation, :translation), 
    Tuple{SVector{3, T}, SMatrix{3, 3, T, 9}, CartesianVector{T}}}
# Transformation rule:
#   1. Scale
#   2. Rotate
#   3. Translate

transform(g::AbstractPrimitive, t::Transformations) =
    transform(transform(transform(g, t.scale), t.rotation), t.translation)

transform(p::AbstractVolumePrimitive, v::CartesianVector) = p + v
transform(p::AbstractVolumePrimitive, r::AbstractMatrix) = r * p
transform(p::AbstractVolumePrimitive, s::SVector{3}) = s * p

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

