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

Edge{T}(a::CylindricalPoint{T}, b::CylindricalPoint{T}) where {T} = Edge{T}(CartesianPoint(a), CartesianPoint(b))

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