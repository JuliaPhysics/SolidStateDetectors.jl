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
    a::CartesianPoint{T}
    b::CartesianPoint{T}
end

function distance_to_line(point::AbstractCoordinatePoint{T}, seg::LineSegment{T})::T where {T}
    p = CartesianPoint(point)

    v_ab = normalize(CartesianVector{T}(seg.b - seg.a))
    v_ap = CartesianVector{T}(p - seg.a)

    proj = dot(v_ab, v_ap)

    if geom_round(proj) ≤ zero(T)
        return norm(seg.a - p)
    end

    v_bp = CartesianVector{T}(p - seg.b)
    if geom_round(dot(v_ab, v_bp)) ≥ zero(T)
        return norm(seg.b - p)
    end

    return sqrt(abs(dot(v_ap, v_ap) - proj^2))
end
