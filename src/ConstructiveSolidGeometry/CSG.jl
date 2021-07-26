"""
    struct CSGUnion{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}
        
a || b
"""
struct CSGUnion{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}
    a::A
    b::B
end

in(p::AbstractCoordinatePoint{T}, csg::CSGUnion; csgtol::T = csg_default_tol(T)) where {T} = 
    in(p, csg.a; csgtol = csgtol) || in(p, csg.b; csgtol = csgtol)
(+)(a::A, b::B) where {T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} = CSGUnion{T,A,B}(a, b)

# read-in
function Geometry(::Type{T}, t::Type{CSGUnion}, v::Vector{<:AbstractDict}, input_units::NamedTuple, transformations::Transformations{T}) where {T}
    sum( broadcast(x-> Geometry(T, x, input_units, transformations), v) )
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

in(p::AbstractCoordinatePoint{T}, csg::CSGIntersection; csgtol::T = csg_default_tol(T)) where {T} = 
    in(p, csg.a; csgtol = csgtol) && in(p, csg.b; csgtol = csgtol)
(&)(a::A, b::B) where {T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} = CSGIntersection{T,A,B}(a, b)

# read-in
function Geometry(::Type{T}, ::Type{CSGIntersection}, v::Vector{<:AbstractDict}, input_units::NamedTuple, transformations::Transformations{T}) where {T}
    parts = broadcast(x-> Geometry(T, x, input_units, transformations), v) 
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

in(p::AbstractCoordinatePoint{T}, csg::CSGDifference; csgtol::T = csg_default_tol(T)) where {T} = 
    in(p, csg.a; csgtol = csgtol) && !in(p, csg.b; csgtol = csgtol)

function (-)(a::A, b::B) where {T, A <: AbstractGeometry{T}, B <: AbstractConstructiveGeometry{T}} 
    ob = switchClosedOpen(b)
    CSGDifference{T,A,typeof(ob)}(a, ob)
end

function (-)(a::A, b::B) where {T, A <: AbstractGeometry{T}, B <: AbstractVolumePrimitive{T}} 
    ob = OpenPrimitive(b)
    CSGDifference{T,A,typeof(ob)}(a, ob)
end

# read-in
function Geometry(::Type{T}, ::Type{CSGDifference}, v::Vector{<:AbstractDict}, input_units::NamedTuple, transformations::Transformations{T}) where {T}
    Geometry(T, v[1], input_units, transformations) - sum(broadcast(x-> Geometry(T, x, input_units, transformations), v[2:end]))
end

Dictionary(g::CSGDifference{T}) where {T} = OrderedDict{String,Any}("difference" => OrderedDict[Dictionary(g.a), Dictionary(g.b)])

function Geometry(::Type{T}, CSG::Type{<:AbstractConstructiveGeometry}, v::Vector{Any}, input_units::NamedTuple, transformations::Transformations{T}) where {T} 
    Geometry(T, CSG, [g for g in v], input_units, transformations)
end

(+)(csg::A, v::CartesianVector) where {A <: AbstractConstructiveGeometry} = A(csg.a + v, csg.b + v)

(*)(r::AbstractMatrix, csg::A) where {A <: AbstractConstructiveGeometry} = A(r * csg.a, r * csg.b)

surfaces(csg::AbstractConstructiveGeometry) = vcat(surfaces(csg.a)..., surfaces(csg.b)...)

function switchClosedOpen(csg::CSGDifference) 
    CSGDifference(switchClosedOpen(csg.a), switchClosedOpen(csg.b))
end
function switchClosedOpen(csg::CSGUnion) 
    CSGUnion(switchClosedOpen(csg.a), switchClosedOpen(csg.b))
end
function switchClosedOpen(csg::CSGIntersection) 
    CSGIntersection(switchClosedOpen(csg.a), switchClosedOpen(csg.b))
end


sample(csg::AbstractConstructiveGeometry{T}, sampling...) where {T} = vcat(sample(csg.a, sampling...), sample(csg.b, sampling...))