"""
    abstract type AbstractLine{T, N, S} <: AbstractGeometry{T, N} end

T: eltype
N: Dimension
S: Coordinate System: `:Cartesian` or `:Cylindrical`
"""
abstract type AbstractLine{T, N, S} <: AbstractGeometry{T, N} end

"""
    struct Line{T, N, S} <: AbstractLine{T, N, S}

<----A----------B---->
"""
struct Line{T, N, S} <: AbstractLine{T, N, S}
    org::AbstractCoordinatePoint{T, N, S}
    dir::CoordinateVector{T, N, S}
end

function Line(p1::CartesianPoint{T}, p2::CartesianPoint{T})::Line{T, 3, :Cartesian} where {T}
    return Line{T, 3, :Cartesian}(p1, p2 - p1)
end
function in(pt::CartesianPoint{T}, l::Line{T, 3, :Cartesian}; atol::Real = T(0))::Bool where {T}
    norm((pt - l.org) × l.dir) <= atol
end

"""
    struct Ray{T, N, S} <: AbstractLine{T, N, S}

    [A----------B---->
"""
struct Ray{T, N, S} <: AbstractLine{T, N, S}
    org::AbstractCoordinatePoint{T, N, S}
    dir::CoordinateVector{T, N, S}
end
function Ray(p1::CartesianPoint{T}, p2::CartesianPoint{T})::Ray{T, 3, :Cartesian} where {T}
    return Ray{T, 3, :Cartesian}(p1, p2 - p1)
end
function in(pt::CartesianPoint{T}, l::Ray{T, 3, :Cartesian}; atol::Real = T(0))::Bool where {T}
    d::CartesianVector{T} = pt - l.org
    if norm(d × l.dir) <= atol
        return d ⋅ l.dir >= 0
    end
    return false
end

"""
    struct LineSegment{T, N, S} <: AbstractLine{T, N, S}

    [A----------B]
"""
struct LineSegment{T, N, S} <: AbstractLine{T, N, S}
    org::AbstractCoordinatePoint{T, N, S}
    dir::CoordinateVector{T, N, S}
end

function LineSegment(p1::CartesianPoint{T}, p2::CartesianPoint{T})::LineSegment{T, 3, :Cartesian} where {T}
    return LineSegment{T, 3, :Cartesian}(p1, p2 - p1)
end
function in(pt::CartesianPoint{T}, l::LineSegment{T, 3, :Cartesian}; atol::T = T(0))::Bool where {T}
    d::CartesianVector{T} = pt - l.org
    if norm(d × l.dir) <= atol
        d2::CartesianVector{T} = d - l.dir
        return (d ⋅ l.dir >= 0) && (d2 ⋅ l.dir <= 0)
    end
    return false
end


function iscolinear(l1::AbstractLine{T, 3, :Cartesian}, l2::AbstractLine{T, 3, :Cartesian}, atol::Real = T(0))::Bool where {T}
    return norm(l1.dir × l2.dir) <= atol # l1.dir == l2.dir || norm(l1.dir × l2.dir) == 0
end



# Plot Recipes:
@recipe function f(l::AbstractLine{T, 3, :Cartesian}) where {T}
    x::Vector{T} = [l.org.x, l.org.x + l.dir.x] 
    y::Vector{T} = [l.org.y, l.org.y + l.dir.y] 
    z::Vector{T} = [l.org.z, l.org.z + l.dir.z] 
    x, y, z
end

@recipe function f(l::AbstractLine{T, 2, :Cartesian}) where {T}
    x::Vector{T} = [l.org.x, l.org.x + l.dir.x] 
    y::Vector{T} = [l.org.y, l.org.y + l.dir.y] 
    x, y
end