
abstract type AbstractContact{T} end

"""
    mutable struct Contact{T, D} <: AbstractContact{T}

T: Type of precision.
D: `:N` for n- and `:P` for p-contacts.
"""
mutable struct Contact{T, D} <: AbstractContact{T}
    potential::T
    id::Int 
    name::String
    geometry::Vector{AbstractGeometry{T}}
end

get_contact_type(c::Contact{T, D}) where {T <: AbstractFloat, D} = D 


@inline function in(pt::StaticVector{3, T}, c::Contact{T})::Bool where {T}
    return in(pt, c.geometry)
end


function Contact{T, D}(dict::Dict{Any, Any}, inputunit::Unitful.Units)::Contact{T, D} where {T <: AbstractFloat, D}
    return Contact{T, D}( dict["potential"], dict["channel"], dict["name"], Geometry(T, dict["geometry"], inputunit ) )
end

function in(p::CylindricalPoint{T}, c::Contact{T}, rs::Vector{T}) where T
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