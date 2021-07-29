abstract type AbstractVirtualVolume{T} end

in(pt::AbstractCoordinatePoint{T, 3}, avv::AbstractVirtualVolume{T}) where {T <: SSDFloat} = in(pt, avv.geometry)


function construct_virtual_volume(T, pass::Dict, input_units::NamedTuple, transformations::Transformations)
    construct_virtual_volume(T, pass, input_units, Val{Symbol(pass["model"])}, transformations)
end

"""
    struct DeadVolume{T, G} <: AbstractVirtualVolume{T}
        
Volume inside which the charge drift is set to zero.

## Parametric types
* `T`: Precision type.
* `G`: Type of `geometry`.

## Fields
* `name::String`: Name of the dead volume, relevant for plotting.
* `geometry::G`: Geometry of the dead volume, see [Constructive Solid Geometry (CSG)](@ref).

## Definition in Configuration File
A `DeadVolume` is defined through an entry in the `virtual_drift_volumes`
array of a detector with `model: dead`. It needs a `geometry` and can optionally be given a `name`.

An example definition of dead volumes looks like this:
```yaml 
virtual_drift_volume:
  - name: Volume 1
    model: dead
    geometry: # ...
  - name: Volume 2
    model: dead
    geometry: # ...
```
"""
struct DeadVolume{T, G} <: AbstractVirtualVolume{T}
    name::String
    geometry::G
end

function modulate_driftvector(sv::CartesianVector{T}, pt::CartesianPoint{T}, tl::DeadVolume{T})::CartesianVector{T} where {T <: SSDFloat}
    return CartesianVector{T}(0,0,0)
end

function DeadVolume{T}(dict::Dict, input_units::NamedTuple, transformations::Transformations) where T <: SSDFloat
    n = haskey(dict, "name") ? dict["name"] : "external part"
    g = Geometry(T, dict["geometry"], input_units, transformations)
    return DeadVolume{T, typeof(g)}(n, g)
end

function construct_virtual_volume(T, pass::Dict, input_units::NamedTuple, ::Type{Val{:dead}}, transformations::Transformations)
    DeadVolume{T}(pass, input_units, transformations)
end



struct ArbitraryDriftModificationVolume{T} <: AbstractVirtualVolume{T}
    name::String
    id::Int
    geometry::AbstractGeometry{T}
end

function modulate_driftvector(sv::CartesianVector{T}, pt::CartesianPoint{T}, tl::ArbitraryDriftModificationVolume{T})::CartesianVector{T} where {T <: SSDFloat}
    modulate_driftvector(sv, pt, tl, Val{tl.id})
end
function modulate_driftvector(sv::CartesianVector{T}, pt::CartesianPoint{T}, tl::ArbitraryDriftModificationVolume{T}, ::Type{Val{id}})::CartesianVector{T} where {T <: SSDFloat, id}
    error("""
        This function needs to be overwritten by the user. Use `::Type{Val{<id>}}` for the last argument. 
        <id> is the corresponding id specified in the configuration file of the detector. 
    """)
end


function ArbitraryDriftModificationVolume{T}(dict::Dict, input_units::NamedTuple, transformations::Transformations) where T <: SSDFloat
    n = haskey(dict, "name") ? dict["name"] : "external part"
    g = Geometry(T, dict["geometry"], input_units, transformations)
    id = Int(dict["id"])
    return ArbitraryDriftModificationVolume{T}(n, id, g)
end

function construct_virtual_volume(T, pass::Dict, input_units::NamedTuple, ::Type{Val{:arbitrary}}, transformations::Transformations)
    ArbitraryDriftModificationVolume{T}(pass, input_units, transformations)
end

