
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

