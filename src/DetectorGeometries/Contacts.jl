
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
    geometry::AbstractGeometry{T}
    geometry_positive::Vector{AbstractGeometry{T}}
    geometry_negative::Vector{AbstractGeometry{T}}
end

get_contact_type(c::Contact{T, D}) where {T <: SSDFloat, D} = D

function Contact{T, D}(dict::Union{Dict{String,Any}, Dict{Any, Any}}, inputunit::Unitful.Units)::Contact{T, D} where {T <: SSDFloat, D}
    haskey(dict, "channel") ? channel = dict["channel"] : channel = -1
    haskey(dict, "material") ? material = material_properties[materials[dict["material"]]] : material = material_properties[materials["HPGe"]]
    geometry, geometry_positive, geometry_negative =  Geometry(T, dict["geometry"], inputunit )
    return Contact{T, D}( dict["potential"], material, channel, dict["name"], geometry, geometry_positive, geometry_negative )
end

@inline function in(pt::StaticVector{3, T}, c::Contact{T})::Bool where {T}
    return in(pt, c.geometry)
end

function in(p::AbstractCoordinatePoint{T}, c::Contact{T}, rs::Vector{T}) where T
    rv = false
    for g in c.geometry_positive
        if typeof(g) in []#[ConeMantle{T}]
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

function find_closest_gridpoint(point::CylindricalPoint{T}, grid::CylindricalGrid{T}) where T
    return Int[searchsortednearest( grid[:r].ticks, point.r), searchsortednearest(grid[:φ].ticks, point.φ), searchsortednearest(grid[:z].ticks, point.z)]
end
function find_closest_gridpoint(point::CartesianPoint{T}, grid::CylindricalGrid{T}) where T
    find_closest_gridpoint(CylindricalPoint(point),grid)
end

function find_closest_gridpoint(point::CartesianPoint{T}, grid::CartesianGrid{T}) where T
    return Int[searchsortednearest( grid[:x].ticks, point.x), searchsortednearest(grid[:y].ticks, point.y), searchsortednearest(grid[:z].ticks, point.z)]
end
function find_closest_gridpoint(point::CylindricalPoint{T}, grid::CartesianGrid{T}) where T
    find_closest_gridpoint(CartesianPoint(point),grid)
end

function paint_contact(contact::Contact, grid::CylindricalGrid{T}) where T
    stepsize::Vector{T}= [minimum(diff(grid[:r].ticks)), IntervalSets.width(grid[:φ].interval) == 0.0 ? 0.5236 : minimum(diff(grid[:φ].ticks)), minimum(diff(grid[:z].ticks))]
    samples = filter(x-> x in contact.geometry, vcat([sample(g, stepsize) for g in contact.geometry_positive]...))
    contact_gridpoints = unique!([find_closest_gridpoint(sample_point,grid) for sample_point in samples])
    return contact_gridpoints
end

function paint_contact(contact::Contact, grid::CartesianGrid{T}) where T
    stepsize::Vector{T}= [minimum(diff(grid[:x].ticks)), minimum(diff(grid[:y].ticks)), minimum(diff(grid[:z].ticks))]
    samples = filter(x-> x in contact.geometry, vcat([sample(g, stepsize) for g in contact.geometry_positive]...))
    contact_gridpoints = unique!([find_closest_gridpoint(sample_point,grid) for sample_point in samples])
    return contact_gridpoints
end
