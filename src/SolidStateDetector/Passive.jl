abstract type AbstractPassive{T} <: AbstractObject{T} end

mutable struct Passive{T} <: AbstractPassive{T}
    name::String
    id::Int
    potential::Union{Symbol,T}
    temperature::Union{T,Missing}
    material::NamedTuple
    geometry::AbstractGeometry{T}
    geometry_positive::Vector{AbstractGeometry{T}}
    geometry_negative::Vector{AbstractGeometry{T}}

    Passive{T}() where T <: SSDFloat = new{T}()
end

function Passive{T}(dict::Dict, inputunit_dict::Dict{String,Unitful.Units}) where T <: SSDFloat
    pass = Passive{T}()
    haskey(dict, "name") ? pass.name = dict["name"] : name = "external part"
    haskey(dict, "id") ? pass.id = dict["id"] : id = -1
    haskey(dict,"potential") ? pass.potential = dict["potential"] : potential = :floating
    haskey(dict, "temperature") ? pass.temperature = dict["temperature"] : pass.temperature = missing
    pass.material = material_properties[materials[dict["material"]]]
    pass.geometry = Geometry(T, dict["geometry"], inputunit_dict)
    pass.geometry_positive, pass.geometry_negative = get_decomposed_volumes(pass.geometry)
    return pass
end
