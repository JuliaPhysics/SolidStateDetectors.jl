abstract type AbstractVirtualVolume{T} end

in(p::AbstractCoordinatePoint{T, 3}, avv::AbstractVirtualVolume{T}) where {T <: SSDFloat} = in(p, avv.geometry)


struct DeadVolume{T} <: AbstractVirtualVolume{T}
    name::String
    geometry::AbstractGeometry{T}
end

function modulate_driftvector(sv::CartesianVector{T}, cp::CartesianPoint{T}, tl::DeadVolume{T})::CartesianVector{T} where {T <: SSDFloat}
    return CartesianVector{T}(0,0,0)
end

function DeadVolume{T}(dict::Dict, input_units::NamedTuple, transformations = missing) where T <: SSDFloat
    n = haskey(dict, "name") ? dict["name"] : "external part"
    g = transform(Geometry(T, dict["geometry"], input_units), transformations)
    return DeadVolume{T}(n, g)
end



struct ArbitraryDriftModificationVolume{T} <: AbstractVirtualVolume{T}
    name::String
    id::Int
    geometry::AbstractGeometry{T}
end

function modulate_driftvector(sv::CartesianVector{T}, cp::CartesianPoint{T}, tl::ArbitraryDriftModificationVolume{T})::CartesianVector{T} where {T <: SSDFloat}
    modulate_driftvector(sv, cp, tl, Val{tl.id})
end
function modulate_driftvector(sv::CartesianVector{T}, cp::CartesianPoint{T}, tl::ArbitraryDriftModificationVolume{T}, ::Type{Val{id}})::CartesianVector{T} where {T <: SSDFloat, id}
    error("""
        This function needs to be overwritten by the user. Use `::Type{Val{<id>}}` for the last argument. 
        <id> is the corresponding id specified in the configuration file of the detector. 
    """)
end


function ArbitraryDriftModificationVolume{T}(dict::Dict, input_units::NamedTuple, transformations = missing) where T <: SSDFloat
    n = haskey(dict, "name") ? dict["name"] : "external part"
    g = transform(Geometry(T, dict["geometry"], input_units), transformations)
    id = Int(dict["id"])
    return ArbitraryDriftModificationVolume{T}(n, id, g)
end

