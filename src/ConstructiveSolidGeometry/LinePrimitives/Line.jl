struct Line{T} <: AbstractLinePrimitive{T}
    origin::CartesianPoint{T}
    direction::CartesianVector{T}
    
    Line{T}(origin, direction) where T = new{T}(origin,direction)
    function Line(origin::PT, direction::DIR) where {PT<:CartesianPoint, DIR<:CartesianVector}
        #Type promotion happens here
        eltypes = _csg_get_promoted_eltype.((PT, DIR))
        T = float(promote_type(eltypes...))
        Line{T}(origin, direction)
    end
end

# function Line(origin::CartesianPoint{T}, direction::CartesianVector{U}) where {T,U}
#     R = promote_type(T,U)
#     return Line(convert(CartesianPoint{R}, origin), convert(CartesianVector{R}, direction))
# end

Line(a::CartesianPoint, b::CartesianPoint) = Line(a, normalize(b - a))

function Line(;
    origin = zero(CartesianPoint{Int}), 
    direction = CartesianVector{Int}(0,0,1)
) 
    Line(origin, direction)
end

function Line{T}(;
    origin = zero(CartesianPoint{Float64}), 
    direction = CartesianVector{Float64}(0,0,1)
) where {T}
    Line{T}(origin, direction)
end

distance(pt::CartesianPoint, l::Line) = norm((pt - l.origin) × l.direction) / norm(l.direction)

function _transform_into_object_coordinate_system(l::Line{T}, p::AbstractPrimitive) where {T}
    origin = _transform_into_object_coordinate_system(l.origin, p) 
    direction = _transform_into_object_coordinate_system(l.direction, p)
    Line( origin, direction )  
end

struct LineSegment{T}
    p1::SVector{3,T}
    p2::SVector{3,T}
end

function distance_to_line(pt::SVector{3,T}, seg::LineSegment{T})::T where T
    v = seg.p2 - seg.p1
    w = pt - seg.p1
    c1 = dot(w,v)
    if c1 <= zero(T)
        return norm(pt - seg.p1)
    end
    c2 = dot(v,v)
    if c2 <= c1
        return norm(pt - seg.p2)
    end
    b = c1 / c2
    closest = seg.p1 + b*v
    return norm(pt - closest)
end
