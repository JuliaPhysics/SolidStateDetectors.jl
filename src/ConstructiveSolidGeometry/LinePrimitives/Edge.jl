struct Edge{T} <: AbstractLinePrimitive{T}
    a::CartesianPoint{T}
    b::CartesianPoint{T}
end

direction(e::Edge) = e.b - e.a

Line(e::Edge{T}) where {T} = Line{T}(e.a, direction(e))

function distance(pt::CartesianPoint{T}, e::Edge{T}) where {T}
    return if (pt - e.b) ⋅ direction(e) >= 0
        norm(pt - e.b)
    elseif (pt - e.a) ⋅ direction(e) <= 0 
        norm(pt - e.a)
    else
        distance(pt, Line(e))
    end
end

function sample(e::Edge{T}; n = 2)::Vector{CartesianPoint{T}} where {T}
    xs = range(zero(T), stop = one(T), length = n)
    pts = Vector{CartesianPoint{T}}(undef, n)
    dir = direction(e)
    for i in eachindex(pts)
        pts[i] = e.a + xs[i] .* dir
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