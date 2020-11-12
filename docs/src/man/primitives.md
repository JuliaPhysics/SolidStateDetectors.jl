# CSG Primitives

CSG = Constructive Solid Geometries

## Primitive: Box

```json
{
    "type": "box",
    "x": {
        "from": 0.0,
        "to":  10.0
    },
    "y":{
        "from": 0.0,
        "to":  10.0
    },
    "z": {
        "from": 0.0,
        "to":  10.0
    }
}
```

## Primitve: Tube

```json
{
    "type": "tube",
    "r": {
        "from": 3.0,
        "to":   7.0
    },
    "phi": {
        "from": 0.0,
        "to": 360.0
    },
    "z": {
        "from": 0.0,
        "to":  10.0
    }
}
```

## Primitive: Cone

```json
{
    "type": "cone",
    "r": {
        "bottom": {
            "from": 10.0,
            "to": 20.0
        },
        "top": {
            "from": 5.0,
            "to": 25.0
        }
    },
    "phi": {
        "from": 0.0,
        "to": 360.0
    },
    "z": {
        "from": 0.0,
        "to":  10.0
    }
}
```

## Primitive: Sphere

```json
{
    "type": "sphere",
    "r": 30,
    "translate": {
        "x": 0,
        "y": 0,
        "z": 0
    }
}
```

## Add new Primitive

If you need a primitive which is not yet implemented you can implemented yourself and make a PR on GitHub :)

There are a few functions you will have to define for a new primitive. Here is a layout for a new primitive.
Place the struct and all the functions in a file `NewPrimitive.jl` inside `<SSD>/src/Geometries/VolumePrimitives`
and include it in the file `<SSD>/src/Geometries/VolumePrimitives/VolumePrimitives.jl`

Note that internally everything is in SI units. So please have a look in the already implemented primitives (`Box`, `Sphere`, `Tube`, ...) how the conversion from config file units to SI units is done. 

```julia
struct NewPrimitive{T} <: AbstractVolumePrimitive{T, 3}
    # fields of the new primitive
end

function in(pt::CartesianPoint{T}, g::NewPrimitive{T})::Bool where {T <: SSDFloat}
    # pt in NewPrimitive?
end
function in(pt::CylindricalPoint{T}, g::NewPrimitive{T})::Bool where {T <: SSDFloat}
    # pt in NewPrimitive?
end

# You also have to implement the function to obtain the primitive from a config file (so an dic)
# You also should provide a example config file containing this new primitive
function NewPrimitive{T}(dict::Union{Dict{Any, Any}, Dict{String, Any}}, inputunit_dict::Dict{String,Unitful.Units})::NewPrimitive{T} where {T <: SSDFloat}
    # ... parse values from dict to NewPrimitive{T}(...)
    return NewPrimitive{T}(
        # ...
    )
end
function Geometry(T::DataType, t::Val{:newprimitive}, dict::Dict{Any, Any}, inputunit_dict::Dict{String,Unitful.Units})
    return NewPrimitive{T}(dict, inputunit_dict)
end

# Also a plot recipe for this new primitive should be provided:
@recipe function f(cb::NewPrimitive{T}) where {T <: SSDFloat}
    label --> "NewPrimitive"
    # ...
end

# For proper grid creation we also need the function get_important_points:
function get_important_points(g::NewPrimitive{T}, ::Val{:r})::Vector{T} where {T <: SSDFloat}
    return T[ # ... ]
end
function get_important_points(g::NewPrimitive{T}, ::Val{:φ})::Vector{T} where {T <: SSDFloat}
    return T[ # ... ]
end
function get_important_points(g::NewPrimitive{T}, ::Val{:z})::Vector{T} where {T <: SSDFloat}
    return T[ # ... ]
end
function get_important_points(g::NewPrimitive{T}, ::Val{:x})::Vector{T} where {T <: SSDFloat}
    return T[ # ... ]
end
function get_important_points(g::NewPrimitive{T}, ::Val{:y})::Vector{T} where {T <: SSDFloat}
    return T[ # ... ]
end

# and a sample function to paint the primitive on the grid (necessary if the object is small)
function sample(cb::NewPrimitive{T}, stepsize::Vector{T})  where {T <: SSDFloat}
    samples  = CartesianPoint{T}[]
    # ...
    return samples
end

# add a (+) method to shift the primitive 
function (+)(b::NewPrimitive{T}, translate::Union{CartesianVector{T},Missing})::NewPrimitive{T} where {T <: SSDFloat}
    # ...
    return NewPrimitive( # ...  )
end
```