const Transformations{T} = NamedTuple{(:rotation, :translation), Tuple{SMatrix{3, 3, T, 9}, CartesianVector{T}}}
# Transformation rule:
#   1. Rotate
#   2. Translate
# Not sure yet if Rotations will always be of type SMatrix{3, 3, T, 9} -> Any

# Transformations{T}() where {T} = (rotation = one(SMatrix{3, 3, T, 9}), translation = zero(CartesianVector{T}))

# rotate(p::P, r::AbstractMatrix) where {P <: AbstractPrimitive} = P(p, origin = r * p.origin, rotation = r * p.rotation)

function rotate(p::P, r::AbstractMatrix) where {P <: AbstractPrimitive}
    # Compose rotation with the primitive's local frame:
    P(p, origin = cartesian_zero + r * (p.origin - cartesian_zero), rotation = r * p.rotation)
end
# ToDo: Make this obsolete and then remove it:
(*)(r::AbstractMatrix, p::AbstractPrimitive) = rotate(p, r)

translate(p::P, v::CartesianVector) where {P <: AbstractPrimitive} = P(p, origin = p.origin + v, rotation = p.rotation)
# ToDo: Make this obsolete and then remove it:
(+)(p::AbstractPrimitive, v::CartesianVector) = translate(p, v)

transform(g::AbstractPrimitive, t::Transformations) =
    translate(rotate(g, t.rotation), t.translation)

transform(csg::G, t::Transformations) where {G <: AbstractConstructiveGeometry} = G(transform(csg.a, t), transform(csg.b, t))


function combine_transformations(current::Transformations, new::Transformations)
    r = new.rotation * current.rotation    
    t = new.rotation * current.translation
    t = t + new.translation
    (rotation = r, translation = t)
end


function Dictionary(m::SMatrix{3,3,T,9}) where {T}
    dict = OrderedDict{String, Any}()
    mat = RotXYZ(m)
    if mat.theta1 == 0 && mat.theta2 == 0
        if mat.theta3 != 0
            dict["Z"] = string(rad2deg(mat.theta3))*"°"
        end
    else
        dict["M"] = m[:]
    end
    return dict
end

Dictionary(pt::CartesianPoint) = [pt.x, pt.y, pt.z]

