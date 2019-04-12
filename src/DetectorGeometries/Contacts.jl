
abstract type AbstractContact{T} end

"""
    mutable struct Contact{T, D} <: AbstractContact{T}

T: Type of precision.
D: `:N` for n- and `:P` for p-contacts.
"""
mutable struct Contact{T, D} <: AbstractContact{T}
    potential::T
    material::NamedTuple
    id::Int
    name::String
    geometry::Vector{AbstractGeometry{T}}
end

get_contact_type(c::Contact{T, D}) where {T <: SSDFloat, D} = D

function Contact{T, D}(dict::Union{Dict{String,Any}, Dict{Any, Any}}, inputunit::Unitful.Units)::Contact{T, D} where {T <: SSDFloat, D}
    haskey(dict, "channel") ? channel = dict["channel"] : channel = -1
    haskey(dict, "material") ? material = material_properties[materials[dict["material"]]] : material = material_properties[materials["HPGe"]]
    return Contact{T, D}( dict["potential"], material, channel, dict["name"], Geometry(T, dict["geometry"], inputunit ) )
end

@inline function in(pt::StaticVector{3, T}, c::Contact{T})::Bool where {T}
    return in(pt, c.geometry)
end

function in(p::AbstractCoordinatePoint{T}, c::Contact{T}, rs::Vector{T}) where T
    rv = false
    for g in c.geometry
        if typeof(g) == ConeMantle{T}
            in(p, g, rs) ? rv = true : nothing
        else
            (p in g) ? rv = true : nothing
        end
    end
    return rv
end

function in(p::AbstractCoordinatePoint{T}, v::AbstractVector{<:AbstractContact{T}}) where T
    rv = false
    for contact in v
        if p in contact
            rv = true
        end
    end
    return rv
end
