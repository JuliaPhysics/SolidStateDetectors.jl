
abstract type AbstractContact{T} <: AbstractObject{T} end

"""
    mutable struct Contact{T} <: AbstractContact{T}

T: Type of precision.
"""
mutable struct Contact{T} <: AbstractContact{T}
    potential::T
    material::NamedTuple
    id::Int
    name::String
    geometry::AbstractGeometry{T}
end


function Contact{T}(dict::Union{Dict{String,Any}, Dict{Any, Any}}, input_units::NamedTuple, transformations::Vector{CSGTransformation})::Contact{T} where {T <: SSDFloat}
    haskey(dict, "channel") ? channel = dict["channel"] : channel = -1
    haskey(dict, "material") ? material = material_properties[materials[dict["material"]]] : material = material_properties[materials["HPGe"]]
    haskey(dict,"name") ? name = dict["name"] : name = ""
    geometry = transform(Geometry(T, dict["geometry"], input_units), transformations)
    return Contact{T}( dict["potential"], material, channel, name, geometry )
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
