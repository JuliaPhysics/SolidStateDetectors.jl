"""
    struct GeometryUnion{T, N, A, B} <: AbstractSet{T, N}

a || b
"""
struct GeometryUnion{T <: Real, N, A, B} <: AbstractSet{T, N}
    # just `Union` seems to conflict with Core.Union... :(
    a::A
    b::B
end

function GeometryUnion{T, N}(a::A, b::B) where {T, N, A, B}
    GeometryUnion{T, N, A, B}(a, b)
end

function in(pt::StaticVector{N, T}, set::GeometryUnion{T, N})::Bool where {T <: Real, N}
    inside::Bool = false
    if in(pt, set.a) || in(pt, set.b)
        inside = true
    end
    return inside
end

"""
    struct Intersection{T, N, A, B} <: AbstractSet{T, N}

a && b
"""
struct Intersection{T, N, A, B} <: AbstractSet{T, N}
    a::A
    b::B
end
function Intersection{T, N}(a::A, b::B) where {T, N, A, B}
    Intersection{T, N, A, B}(a, b)
end
function in(pt::StaticVector{N, T}, set::Intersection{T, N})::Bool where {T <: Real, N}
    inside::Bool = false
    if in(pt, set.a) && in(pt, set.b)
        inside = true
    end
    return inside
end

"""
    struct Difference{T, N, A, B} <: AbstractSet{T, N}

a && !b
"""
struct Difference{T, N, A, B} <: AbstractSet{T, N}
    a::A
    b::B
end
function Difference{T, N}(a::A, b::B) where {T, N, A, B}
    Difference{T, N, A, B}(a, b)
end
function in(pt::StaticVector{N, T}, set::Difference{T, N})::Bool where {T <: Real, N}
    inside::Bool = false
    if in(pt, set.a) && !in(pt, set.b)
        inside = true
    end
    return inside
end


function (+)(a::A, b::B)::AbstractSet{T, N} where {T, N, A <: AbstractGeometry{T, N}, B <: AbstractGeometry{T, N}}
    return GeometryUnion{T, N}(a, b)
end

function (-)(a::A, b::B)::AbstractSet{T, N} where {T, N, A <: AbstractGeometry{T, N}, B <: AbstractGeometry{T, N}}
    return Difference{T, N}(a, b)
end

function (&)(a::A, b::B)::AbstractSet{T, N} where {T, N, A <: AbstractGeometry{T, N}, B <: AbstractGeometry{T, N}}
    return Intersection{T, N}(a, b)
end
