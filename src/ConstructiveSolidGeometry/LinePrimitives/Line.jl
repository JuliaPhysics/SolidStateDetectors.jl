struct Line{T} <: AbstractLinePrimitive{T}
    origin::CartesianPoint{T}
    direction::CartesianVector{T}
end

function distance(A::CartesianPoint, l::Line)
    B = l.origin 
    C = B + l.direction
    return norm((A - B) × l.direction) / norm(l.direction)
end
