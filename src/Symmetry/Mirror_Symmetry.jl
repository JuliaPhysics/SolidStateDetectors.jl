struct MirrorSymmetry{T}
    symmetry_plane::Plane{T} #Plane{T}(origin, normal)
end

function MirrorSymmetry{T}(origin::CylindricalPoint{T}, normal::CylindricalVector{T}) where {T}
    _origin = CartesianPoint(origin)
    _normal = CartesianVector{T}(cos(origin[2]) * normal[1] - sin(origin[2]) * normal[2], 
        sin(origin[2]) * normal[1] + cos(origin[2]) * normal[2], normal[3])
    MirrorSymmetry{T}(Plane{T}(_origin, _normal))
end

function MirrorSymmetry{T}(origin::CartesianPoint{T}, normal::CartesianVector{T}) where {T}
    MirrorSymmetry{T}(Plane{T}(origin, normal))
end

function MirrorSymmetry(axis, value::T, units::NamedTuple = default_unit_tuple()) where {T<:SSDFloat}
    length_unit = units.length
    angle_unit = units.angle
    if axis == "phi" || axis == "φ"
        MirrorSymmetry{T}(CylindricalPoint{T}(r=to_internal_units(1 * length_unit), φ = to_internal_units(value * angle_unit)),
            CylindricalVector{T}(0, to_internal_units(1 * length_unit), 0))
    elseif axis == "x"
        MirrorSymmetry{T}(CartesianPoint{T}(x = to_internal_units(value * length_unit)), CartesianVector{T}(1,0,0))
    elseif axis == "y"
        MirrorSymmetry{T}(CartesianPoint{T}(y = to_internal_units(value * length_unit)), CartesianVector{T}(0,1,0))
    elseif axis == "z"
        MirrorSymmetry{T}(CartesianPoint{T}(z = to_internal_units(value * length_unit)), CartesianVector{T}(0,0,1))
    elseif axis == "r"
        @error "Mirror symmetry along r not defined"
    else
        @error "Wrong axis"
    end
end
        
        