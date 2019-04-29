# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).
"""
    struct RectangularCuboid{T} <: AbstractGeometry{T, 3}

...
"""
struct RectangularCuboid{T} <: AbstractGeometry{T, 3}
    org::CartesianPoint{T}
    e1::CartesianVector{T}
    e2::CartesianVector{T}
    e3::CartesianVector{T}
end

function in(pt::CartesianPoint{T}, rc::RectangularCuboid{T})::Bool where {T <: Real}
    shift::CartesianVector{T} = pt - rc.org
    return  shift ⋅ rc.e1 >= 0 &&
            shift ⋅ rc.e2 >= 0 &&
            shift ⋅ rc.e3 >= 0 &&
            (shift - rc.e1) ⋅ rc.e1 <= 0 &&
            (shift - rc.e2) ⋅ rc.e2 <= 0 &&
            (shift - rc.e3) ⋅ rc.e3 <= 0 
end

function vertices(rc::RectangularCuboid{T})::Vector{CartesianPoint{T}} where {T}
    v::Vector{CartesianPoint{T}} = CartesianPoint{T}[
        rc.org,
        rc.org + rc.e1,
        rc.org + rc.e1 + rc.e2,
        rc.org + rc.e2,
        rc.org + rc.e3,
        rc.org + rc.e3 + rc.e1,
        rc.org + rc.e3 + rc.e1 + rc.e2,
        rc.org + rc.e3 + rc.e2,
    ]    
    return v
end

function LineSegments(rc::RectangularCuboid{T})::Vector{LineSegment{T, 3, :cartesian}} where {T}
    return LineSegment{T, 3, :cartesian}[
        LineSegment{T, 3, :cartesian}( rc.org, rc.e1),
        LineSegment{T, 3, :cartesian}( rc.org + rc.e1,  rc.e2),
        LineSegment{T, 3, :cartesian}( rc.org + rc.e1 + rc.e2, -rc.e1),
        LineSegment{T, 3, :cartesian}( rc.org + rc.e2, -rc.e2),
        LineSegment{T, 3, :cartesian}( rc.org, rc.e3),
        LineSegment{T, 3, :cartesian}( rc.org + rc.e1,  rc.e3),
        LineSegment{T, 3, :cartesian}( rc.org + rc.e1 + rc.e2, rc.e3),
        LineSegment{T, 3, :cartesian}( rc.org + rc.e2, rc.e3),
        LineSegment{T, 3, :cartesian}( rc.org + rc.e3, rc.e1),
        LineSegment{T, 3, :cartesian}( rc.org + rc.e3 + rc.e1,  rc.e2),
        LineSegment{T, 3, :cartesian}( rc.org + rc.e3 + rc.e1 + rc.e2, - rc.e1),
        LineSegment{T, 3, :cartesian}( rc.org + rc.e3 + rc.e2, -rc.e2)
    ]
end