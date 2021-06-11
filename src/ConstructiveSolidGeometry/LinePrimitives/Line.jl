struct Line{T} <: AbstractLinePrimitive{T}
    origin::CartesianPoint{T}
    normal::CartesianVector{T}
end

function distance(A::CartesianPoint, l::Line)
    B = l.origin 
    C = B + l.normal
    return norm((A - B) × l.normal) / norm(l.normal)
end
