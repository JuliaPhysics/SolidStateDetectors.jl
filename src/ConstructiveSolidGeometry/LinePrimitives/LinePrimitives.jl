struct LineSegments{T} <: AbstractLinePrimitive{T}
    points::Vector{CartesianPoint{T}}
end

function (*)(pl::LineSegments{T}, s::SVector{3,T})::LineSegments{T} where {T}
    LineSegments(map(p -> p .* s, pl.points))
end
(*)(s::SVector{3,T}, pl::LineSegments{T}) where {T} = pl * s

function (+)(pl::LineSegments{T}, t::CartesianVector{T})::LineSegments{T} where {T}
    LineSegments(map(p -> p + t, pl.points))
end

function (*)(r::RotMatrix{3,T,TT}, pl::LineSegments{T})::LineSegments{T} where {T, TT}
    LineSegments(map(p -> r * p, pl.points))
end