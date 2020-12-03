struct LineSegments{T} <: AbstractLinePrimitive{T}
    points::Vector{CartesianPoint{T}}
end
