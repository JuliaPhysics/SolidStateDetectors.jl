struct Edge{T} <: AbstractLinePrimitive{T}
    a::CartesianPoint{T}
    b::CartesianPoint{T}
    
    Edge{T}(a, b) where T = new{T}(a, b)
    function Edge(a::A, b::B) where {A<:CartesianPoint, B<:CartesianPoint}
        #Type promotion happens here
        eltypes = _csg_get_promoted_eltype.((A,B))
        T = float(promote_type(eltypes...))
        new{T}(a,b)
    end
end

function Edge(;
    a = CartesianPoint{Int}(0,0,0), 
    b = CartesianPoint{Int}(0,0,1)
) 
    Edge(a, b)
end

function Edge{T}(;
    a = CartesianPoint{Int}(0,0,0), 
    b = CartesianPoint{Int}(0,0,1)
) where {T}
    Edge{T}(a, b)
end

Line(e::Edge) = Line(e.a, e.b)

function distance(pt::CartesianPoint{T}, e::Edge{T}) where {T}
    v = e.b - e.a
    return if (pt - e.b) ⋅ v >= 0
        norm(pt - e.b)
    elseif (pt - e.a) ⋅ v <= 0 
        norm(pt - e.a)
    else
        distance(pt, Line(e))
    end
end

function sample(e::Edge{T}; n = 2)::Vector{CartesianPoint{T}} where {T}
    xs = range(zero(T), stop = one(T), length = n)
    pts = Vector{CartesianPoint{T}}(undef, n)
    v = e.b - e.a
    for i in eachindex(pts)
        pts[i] = e.a + xs[i] .* v
    end
    pts
end
function sample(es::Vector{Edge{T}}; n = 2)::Vector{CartesianPoint{T}} where {T}
    pts = Vector{CartesianPoint{T}}(undef, length(es) * n)
    for i in eachindex(es)
        pts[(i-1)*n+1:i*n] = sample(es[i], n = n)
    end
    pts
end


function distance_to_line(point::AbstractCoordinatePoint{T}, edge::Edge{T})::T where {T}
    return distance(CartesianPoint(point), edge)
    #=
    # ToDo: how is distance_to_line different from distance ? Can this be combined ?
    point = CartesianPoint(point)
    v12 = normalize(CartesianVector{T}(edge.b - edge.a))
    v_point_1 = CartesianVector{T}(point - edge.a)
    proj_on_v12 = dot(v12,v_point_1)
    if geom_round(proj_on_v12) ≤ T(0)
        return norm(edge.a - point)
    else
        v_point_2 = CartesianVector{T}(point - edge.b)
        if geom_round(dot(v12,v_point_2)) ≥ T(0)
            return norm(edge.b - point)
        else
            return sqrt(abs(dot(v_point_1,v_point_1) - proj_on_v12^2))
        end
    end
    =#
end