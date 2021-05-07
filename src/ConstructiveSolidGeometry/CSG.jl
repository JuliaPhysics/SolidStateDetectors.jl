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

Geometry(::Type{T}, CSG::Type{<:AbstractConstructiveGeometry}, v::Vector{Any}, input_units::NamedTuple) where {T} = Geometry(T, CSG, [g for g in v], input_units)
