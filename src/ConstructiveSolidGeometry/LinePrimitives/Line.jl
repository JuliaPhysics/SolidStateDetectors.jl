struct Line{T} <: AbstractLinePrimitive{T}
    origin::CartesianPoint{T}
    direction::CartesianVector{T}
end

distance(pt::CartesianPoint, l::Line) = norm((pt - l.origin) × l.direction) / norm(l.direction)

function _transform_into_object_coordinate_system(l::Line{T}, p::AbstractPrimitive) where {T}
    origin = _transform_into_object_coordinate_system(l.origin, p) 
    direction = rotation(p) * l.direction
    Line( origin, CartesianVector(direction) )  
end