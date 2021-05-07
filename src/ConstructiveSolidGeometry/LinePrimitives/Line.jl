struct Line{T,TP,TL} <: AbstractLinePrimitive{T}
    p1::TP
    p2::TP
    linetype::TL
    function Line( ::Type{T},
                   p1::PT,
                   p2::PT,
                   linetype::Union{Val{:inf}, Val{:ray}, Val{:seg}}
                 ) where {T, PT <: Union{CartesianPoint{T}, PlanarPoint{T}}}
        new{T,typeof(p1),typeof(linetype)}(p1, p2)
    end
end

function LineSegment( ::Type{T},
               p1::PT,
               p2::PT) where {T, PT <: Union{CartesianPoint{T}, PlanarPoint{T}}} #if LineSegment endpoints are p1, and p2
    Line(T, p1, p2, Val(:seg))
end

get_line_vector(l::Line{T,PlanarPoint{T}}) where {T} = normalize(PlanarVector{T}(l.p2 - l.p1))
get_line_vector(l::Line{T,CartesianPoint{T}}) where {T} = normalize(CartesianVector{T}(l.p2 - l.p1))

function _distance_to_line(point::TP, l::Line{T, TP}) where {T, TP <: Union{PlanarPoint{T}, CartesianPoint{T}}}
    v12 = get_line_vector(l)
    v_point_1 = point - l.p1
    proj_on_v12 = dot(v12,v_point_1)
    return v12, v_point_1, proj_on_v12
end

function distance_to_line(point::TP, l::Line{T, TP, Val{:inf}})::T where {T, TP <: Union{PlanarPoint{T}, CartesianPoint{T}}}
    v12, v_point_1, proj_on_v12 = _distance_to_line(point, l)
    return sqrt(abs(dot(v_point_1,v_point_1) - proj_on_v12^2))
end

function distance_to_line(point::TP, l::Line{T, TP, Val{:ray}})::T where {T, TP <: Union{PlanarPoint{T}, CartesianPoint{T}}}
    v12, v_point_1, proj_on_v12 = _distance_to_line(point, l)
    if proj_on_v12 ≤ T(0)
        return norm(l.p1 - point)
    else
        v_point_2 = point - l.p2
        if dot(v12,v_point_2) ≥ T(0)
            return norm(l.p2 - point)
        else
            return sqrt(abs(dot(v_point_1,v_point_1) - proj_on_v12^2))
        end
    end
end

function distance_to_line(point::TP, l::Line{T, TP, Val{:seg}})::T where {T, TP <: Union{PlanarPoint{T}, CartesianPoint{T}}}
    v12, v_point_1, proj_on_v12 = _distance_to_line(point, l)
    if proj_on_v12 ≤ T(0)
        return norm(l.p1 - point)
    else
        v_point_2 = point - l.p2
        if dot(v12,v_point_2) ≥ T(0)
            return norm(l.p2 - point)
        else
            return sqrt(abs(dot(v_point_1,v_point_1) - proj_on_v12^2))
        end
    end
end
