"""
    struct CSGUnion{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}
        
a || b
"""
struct CSGUnion{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}
    a::A
    b::B
end

in(p::AbstractCoordinatePoint, csg::CSGUnion) = in(p, csg.a) || in(p, csg.b)

"""
    struct CSGIntersection{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}

a && b
"""
struct CSGIntersection{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}
    a::A
    b::B
end

in(p::AbstractCoordinatePoint, csg::CSGIntersection) = in(p, csg.a) && in(p, csg.b)


"""
    struct CSGDifference{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}

a && !b
"""
struct CSGDifference{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}
    a::A
    b::B
end

in(p::AbstractCoordinatePoint, csg::CSGDifference) = in(p, csg.a) && !in(p, csg.b)

