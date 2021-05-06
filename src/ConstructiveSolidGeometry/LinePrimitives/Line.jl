struct Line{T,TP,TL} <: AbstractLinePrimitive{T}
    p1::TP
    p2::TP
    linetype::TL
    function Line( ::Type{T},
                   p1::Union{CartesianPoint{T}, PlanarPoint{T}},
                   p2::Union{CartesianPoint{T}, PlanarPoint{T}},
                   linetype::Union{Val{:inf}, Val{:ray}, Val{:seg}}) where {T}
        new{T,typeof(p1),typeof(linetype)}(p1, p2)
    end
end

function Line( ::Type{T},
               p1::Union{CartesianPoint{T}, PlanarPoint{T}},
               p2::Union{CartesianPoint{T}, PlanarPoint{T}}) where {T} #if (inf)Line there are no endpoints
    Line(T, p1, p2, Val(:inf))
end

function Ray( ::Type{T},
               p1::Union{CartesianPoint{T}, PlanarPoint{T}},
               p2::Union{CartesianPoint{T}, PlanarPoint{T}}) where {T} #if Ray endpoint is p1
    Line(T, p1, p2, Val(:ray))
end

function LineSegment( ::Type{T},
               p1::Union{CartesianPoint{T}, PlanarPoint{T}},
               p2::Union{CartesianPoint{T}, PlanarPoint{T}}) where {T} #if LineSegment endpoints are p1, and p2
    Line(T, p1, p2, Val(:seg))
end

Line(l::Line{T}) where {T} = Line(T, l.p1, l.p2)

get_line_vector(l::Line{T,PlanarPoint{T}}) where {T} = normalize(PlanarVector{T}(l.p2 - l.p1))
get_line_vector(l::Line{T,CartesianPoint{T}}) where {T} = normalize(CartesianVector{T}(l.p2 - l.p1))


function distance_to_line(point::Union{PlanarPoint{T}, CartesianPoint{T}}, l::Line{T, <:Any, Val{:inf}})::T where {T}
    v12 = normalize(l.p2 - l.p1)
    v_point_1 = point - l.p1
    proj_on_v12 = dot(v12,v_point_1)
    return sqrt(abs(dot(v_point_1,v_point_1) - proj_on_v12^2))
end

function distance_to_line(point::Union{PlanarPoint{T}, CartesianPoint{T}}, l::Line{T, <:Any, Val{:ray}})::T where {T}
    v12 = get_line_vector(l)
    v_point_1 = point - l.p1
    proj_on_v12 = dot(v12,v_point_1)
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

function distance_to_line(point::Union{PlanarPoint{T}, CartesianPoint{T}}, l::Line{T, <:Any, Val{:seg}})::T where {T}
    v12 = get_line_vector(l)
    v_point_1 = point - l.p1
    proj_on_v12 = dot(v12,v_point_1)
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
