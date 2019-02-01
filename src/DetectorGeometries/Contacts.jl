
abstract type AbstractContact{T} end

"""
    mutable struct Contact{T, D} <: AbstractContact{T}

T: Type of precision.
D: `:N` for n- and `:P` for p-contacts.
"""
mutable struct Contact{T, D} <: AbstractContact{T}
    potential::T
    id::Int 
    geometry::AbstractGeometry{T}
end

@inline function in(pt::StaticVector{3, T}, c::Contact{T})::Bool where {T}
    return in(pt, c.geometry)
end


function Contact{T, D}(dict::Dict{String, Any}, inputunit::Unitful.Units)::Contact{T, D} where {T <: AbstractFloat, D}
    return Contact{T, D}( dict["potential"], dict["id"], CartesianBox3D{T}( dict, inputunit ) )
end