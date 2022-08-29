struct CartesianMirrorSymmetry{T}
    symmetry_plane::Plane{T} #Plane{T}(origin, normal)
end

struct CylindricalMirrorSymmetry{T}
    origin::CylindricalPoint{T}
    normal::CylindricalVector{T}
end

function CartesianMirrorSymmetry{T}(origin::CartesianPoint{T}, normal::CartesianVector{T}) where {T}
    CartesianMirrorSymmetry{T}(Plane{T}(origin, normal))
end

function CartesianMirrorSymmetry(axis, value::T, units::NamedTuple = default_unit_tuple()) where {T<:SSDFloat}
    length_unit = units.length
    angle_unit = units.angle
    if axis == "x"
        CartesianMirrorSymmetry{T}(CartesianPoint{T}(x = to_internal_units(value * length_unit)), CartesianVector{T}(1,0,0))
    elseif axis == "y"
        CartesianMirrorSymmetry{T}(CartesianPoint{T}(y = to_internal_units(value * length_unit)), CartesianVector{T}(0,1,0))
    elseif axis == "z"
        CartesianMirrorSymmetry{T}(CartesianPoint{T}(z = to_internal_units(value * length_unit)), CartesianVector{T}(0,0,1))
    else
        @error "Wrong axis"
    end
end

function CylindricalMirrorSymmetry(axis, value::T, units::NamedTuple = default_unit_tuple()) where {T<:SSDFloat}
    length_unit = units.length
    angle_unit = units.angle
    if axis == "phi" || axis == "φ"
        CylindricalMirrorSymmetry{T}(CylindricalPoint{T}(r=to_internal_units(1 * length_unit), φ = to_internal_units(value * angle_unit)),
            CylindricalVector{T}(0, to_internal_units(1 * length_unit), 0))   
    elseif axis == "z"
        CylindricalMirrorSymmetry{T}(CartesianPoint{T}(z = to_internal_units(value * length_unit)), CartesianVector{T}(0,0,1))
    elseif axis == "r"
        @error "Mirror symmetry along r not defined"
    else
        @error "Wrong axis"
    end 
end
    