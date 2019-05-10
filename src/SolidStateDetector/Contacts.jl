
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
    geometry_positive::Vector{AbstractGeometry{T}}
    geometry_negative::Vector{AbstractGeometry{T}}
end


function Contact{T}(dict::Union{Dict{String,Any}, Dict{Any, Any}}, inputunit_dict::Dict{String,Unitful.Units})::Contact{T} where {T <: SSDFloat}
    haskey(dict, "channel") ? channel = dict["channel"] : channel = -1
    haskey(dict, "material") ? material = material_properties[materials[dict["material"]]] : material = material_properties[materials["HPGe"]]
    haskey(dict,"name") ? name = dict["name"] : name = ""
    geometry =  Geometry(T, dict["geometry"], inputunit_dict )
    geometry_positive, geometry_negative = get_decomposed_volumes(geometry)
    return Contact{T}( dict["potential"], material, channel, name, geometry, geometry_positive, geometry_negative )
end

function println(io::IO, d::Contact{T}) where {T <: SSDFloat}
    println("\t________"*"Contact $(d.id)"*"________\n")
    println("\t---General Properties---")
    println("\t-Potential: \t\t $(d.potential)")
    println("\t-Contact Material: \t $(d.material.name)")
    println()

end


function show(io::IO, d::Contact{T}) where {T <: SSDFloat} println(d) end
function print(io::IO, d::Contact{T}) where {T <: SSDFloat} println(d) end
function display(io::IO, d::Contact{T} ) where {T <: SSDFloat} println(d) end
function show(io::IO,::MIME"text/plain", d::Contact) where {T <: SSDFloat}
    show(io, d)
end
