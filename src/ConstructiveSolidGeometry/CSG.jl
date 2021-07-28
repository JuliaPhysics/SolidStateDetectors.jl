"""
    struct CSGUnion{T, A <: AbstractGeometry{T}, B <: AbstractGeometry{T}} <: AbstractConstructiveGeometry{T}
        
A `CSGUnion` of two geometries `a` and `b` is defined as the set of points that are in at least 
one of either `a` or `b` (`a || b`).

![CSGUnion](../../docs/src/assets/CSGUnion.png)

## Parametric types
* `T`: Precision type.
* `A`: Type of geometry `a`.
* `B`: Type of geometry `b`.

## Fields 
* `a::A`: First geometry to build the union.
* `b::B`: Second geometry to build the union.


## Definition in Configuration File

A `CSGUnion` is defined in the configuration file as part of the `geometry` field 
of an object through the field `union`, followed by an array of geometries from
which the union is constructed. 

An example definition of a `CSGUnion` looks like this:
```yaml
union: # a || b
  - tube: # a
      r: 2
      h: 1
  - tube: # b
      r: 1
      h: 1.5
      origin: 
        z: 0.5
```

!!! note
    If more than two geometries are passed, the `union` is constructed from all of them.
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
        
A `CSGIntersection` of two geometries `a` and `b` is defined as the set of points 
that are both in `a` and in `b` (`a && b`).

![CSGIntersection](../../docs/src/assets/CSGIntersection.png)

## Parametric types
* `T`: Precision type.
* `A`: Type of geometry `a`.
* `B`: Type of geometry `b`.

## Fields 
* `a::A`: First geometry to build the intersection.
* `b::B`: Second geometry to build the intersection.


## Definition in Configuration File

A `CSGIntersection` is defined in the configuration file as part of the `geometry` field 
of an object through the field `intersection`, followed by an array of geometries from
which the intersection is constructed. 

An example definition of a `CSGIntersection` looks like this:
```yaml
intersection: # a && b
  - tube: # a
      r: 2
      h: 1
  - tube: # b
      r: 1
      h: 1.5
      origin: 
        z: 0.5
```

!!! note
    If more than two geometries are passed, the `intersection` is constructed from all of them.
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