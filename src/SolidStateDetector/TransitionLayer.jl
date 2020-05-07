struct DeadVolume{T} <: AbstractVirtualVolume{T}
    name::String
    geometry::AbstractGeometry{T}
end

function modulate_driftvector(sv::CartesianVector{T}, cp::CartesianPoint{T}, tl::DeadVolume{T})::CartesianVector{T} where {T <: SSDFloat}
    return CartesianVector{T}(0,0,0)
end

function DeadVolume{T}(dict::Dict, inputunit_dict::Dict{String,Unitful.Units}) where T <: SSDFloat
    n = haskey(dict, "name") ? dict["name"] : "external part"
    g = Geometry(T, dict["geometry"], inputunit_dict)
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
        <id> is the corresponding id specified in the json-config file of the detector. 
    """)
end


function ArbitraryDriftModificationVolume{T}(dict::Dict, inputunit_dict::Dict{String,Unitful.Units}) where T <: SSDFloat
    n = haskey(dict, "name") ? dict["name"] : "external part"
    g = Geometry(T, dict["geometry"], inputunit_dict)
    id = Int(dict["id"])
    return ArbitraryDriftModificationVolume{T}(n, id, g)
end

