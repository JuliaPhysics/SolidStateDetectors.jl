"""
    struct CSGUnion{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}
        
a || b
"""
struct CSGUnion{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}
    a::A
    b::B
end

in(p::AbstractCoordinatePoint, csg::CSGUnion) = in(p, csg.a) || in(p, csg.b)
(+)(a::A, b::B) where {T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} = CSGUnion{T,A,B}(a, b)

# read-in
function Geometry(::Type{T}, t::Type{CSGUnion}, v::Vector{<:AbstractDict}, input_units::NamedTuple) where {T}
    sum( broadcast(x-> Geometry(T, x, input_units), v) )
end

Dictionary(g::CSGUnion{T}) where {T} = OrderedDict{String,Any}("union" => vcat(UnionDictionary(g.a), UnionDictionary(g.b)))
UnionDictionary(g::CSGUnion{T}) where {T} = vcat(UnionDictionary(g.a), UnionDictionary(g.b))
UnionDictionary(g::AbstractGeometry{T}) where {T} = OrderedDict[Dictionary(g)]

"""
    struct CSGIntersection{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}

a && b
"""
struct CSGIntersection{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}
    a::A
    b::B
end

in(p::AbstractCoordinatePoint, csg::CSGIntersection) = in(p, csg.a) && in(p, csg.b)
(&)(a::A, b::B) where {T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} = CSGIntersection{T,A,B}(a, b)

# read-in
function Geometry(::Type{T}, ::Type{CSGIntersection}, v::Vector{<:AbstractDict}, input_units::NamedTuple) where {T}
    parts = broadcast(x-> Geometry(T, x, input_units), v) 
    reduce(&, parts)
end

Dictionary(g::CSGIntersection{T}) where {T} = OrderedDict{String,Any}("intersection" => vcat(IntersectionDictionary(g.a), IntersectionDictionary(g.b)))
IntersectionDictionary(g::CSGIntersection{T}) where {T} = vcat(IntersectionDictionary(g.a), IntersectionDictionary(g.b))
IntersectionDictionary(g::AbstractGeometry{T}) where {T} = [Dictionary(g)]

"""
    struct CSGDifference{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}

a && !b
"""
struct CSGDifference{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}
    a::A
    b::B
end

in(p::AbstractCoordinatePoint, csg::CSGDifference) = in(p, csg.a) && !in(p, csg.b)
(-)(a::A, b::B) where {T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} = CSGDifference{T,A,B}(a, b)

# read-in
function Geometry(::Type{T}, ::Type{CSGDifference}, v::Vector{<:AbstractDict}, input_units::NamedTuple) where {T}
    Geometry(T, v[1], input_units) - sum( broadcast(x-> Geometry(T, x, input_units), v[2:end]) )
end

Dictionary(g::CSGDifference{T}) where {T} = OrderedDict{String,Any}("difference" => OrderedDict[Dictionary(g.a), Dictionary(g.b)])

Geometry(::Type{T}, CSG::Type{<:AbstractConstructiveGeometry}, v::Vector{Any}, input_units::NamedTuple) where {T} = Geometry(T, CSG, [g for g in v], input_units)

(+)(csg::A, v::CartesianVector) where {A <: AbstractConstructiveGeometry} = A(csg.a + v, csg.b + v)

(*)(r::AbstractMatrix, csg::A) where {A <: AbstractConstructiveGeometry} = A(r * csg.a, r * csg.b)

scale(csg::A, s::SVector{3}) where {A <: AbstractConstructiveGeometry} = A(scale(csg.a, s), scale(csg.b, s))
(*)(s::SVector{3}, csg::A) where {A <: AbstractConstructiveGeometry} = A(scale(csg.a, s), scale(csg.b, s))
(*)(csg::A, s::SVector{3}) where {A <: AbstractConstructiveGeometry} = A(scale(csg.a, s), scale(csg.b, s))

faces(csg::AbstractConstructiveGeometry) = vcat(faces(csg.a)..., faces(csg.b)...)