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
function Geometry(T::DataType, t::Val{:union}, dict::Union{Dict{String,Any}, Dict{Any,Any}}, inputunit_dict::Dict{String,Unitful.Units})
    sum( map(x-> Geometry(T, x, inputunit_dict), dict["parts"]) )
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
(&)(a::A, b::B) where {T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} = Intersection{T,A,B}(a, b)

# read-in
function Geometry(T::DataType, t::Val{:intersection}, dict::Union{Dict{String,Any}, Dict{Any,Any}}, inputunit_dict::Dict{String,Unitful.Units})
    parts = map(x-> Geometry(T, x, inputunit_dict), dict["parts"]) 
    intersection = parts[1]
    for part in parts[2:end] intersection &= part end
    intersection
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
function Geometry(T::DataType, t::Val{:difference}, dict::Union{Dict{String,Any}, Dict{Any,Any}}, inputunit_dict::Dict{String,Unitful.Units})
    Geometry(T, dict["parts"][1], inputunit_dict) - sum( map(x-> Geometry(T, x, inputunit_dict), dict["parts"][2:end]) )
end

