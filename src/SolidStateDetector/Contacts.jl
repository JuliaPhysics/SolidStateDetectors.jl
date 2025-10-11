
abstract type AbstractContact{T} <: AbstractObject{T} end

"""
    struct Contact{T, G, MT} <: AbstractContact{T}
        
Contact of a [`SolidStateDetector`](@ref).

For the simulation of the [`ElectricPotential`](@ref), all contacts are fixed to a constant potential value.

## Parametric types:
* `T`: Precision type.
* `G`: Type of `geometry`.
* `MT`: Type of `material`.

## Fields
* `potential::T`: Potential (in V) to which the contact will be fixed during the calculation of the [`ElectricPotential`](@ref).
* `material::MT`: Material of the contact.
* `id::Int`: Unique id that will unambiguously identify the contact.
* `name::String`: Custom name for the contact, relevant for plotting.
* `geometry::G`: Geometry of the contact, see [Constructive Solid Geometry (CSG)](@ref).

## Definition in Configuration File

A `Contact` is defined in the configuration file through an entry in the `contacts` array of a detector.
It needs `id`, `potential` and `geometry` and can optionally be given a `name` and `material`.

An example definition of contacts looks like this:
```yaml 
contacts:
  - name: "n+ contact"
    id: 1
    potential: 5000V
    material: HPGe # optional
    geometry: # ....
  - name: "p+ contact"
    id: 2
    potential: 0
    material: HPGe #optional
    geometry: # ....
```
"""
struct Contact{T,G,MT} <: AbstractContact{T}
    potential::T
    material::MT
    id::Int
    name::String
    geometry::G
end


function Contact{T}(dict::AbstractDict, input_units::NamedTuple, outer_transformations)::Contact{T} where {T <: SSDFloat}
    id::Int = haskey(dict, "id") ? dict["id"] : -1
    material = haskey(dict, "material") ? material_properties[materials[dict["material"]]] : material_properties[materials["HPGe"]]
    name = haskey(dict,"name") ? dict["name"] : ""
    inner_transformations = parse_CSG_transformation(T, dict, input_units)
    transformations = combine_transformations(inner_transformations, outer_transformations)
    geometry = Geometry(T, dict["geometry"], input_units, transformations)
    return Contact( _parse_value(T, dict["potential"], input_units.potential), material, id, name, geometry )
end

function println(io::IO, d::Contact) 
    println("\t________"*"Contact $(d.id)"*"________\n")
    println("\t---General Properties---")
    println("\t-Potential: \t\t $(d.potential) V")
    println("\t-Contact Material: \t $(d.material.name)")
    println()
end
print(io::IO, d::Contact{T}) where {T} = print(io, "Contact $(d.id) - $(d.potential) V - $(d.material.name)")
show(io::IO, d::Contact) = print(io, d)
show(io::IO,::MIME"text/plain", d::Contact) = show(io, d)
