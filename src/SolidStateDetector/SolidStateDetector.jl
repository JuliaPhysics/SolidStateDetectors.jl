"""
    mutable struct SolidStateDetector{T <: SSDFloat} <: AbstractConfig{T}

"""
struct SolidStateDetector{T,SC,CT,PT,VDM} <: AbstractConfig{T}
    name::String  # optional
    semiconductor::SC
    contacts::CT
    passives::PT
    virtual_drift_volumes::VDM
    
    SolidStateDetector{T}(n::AbstractString,s::SC,c::C,p::P,v::VDM) where {T,SC,C,P,VDM}= new{T,SC,C,P,VDM}(n,s,c,p,v)
end


get_precision_type(::SolidStateDetector{T}) where {T} = T

function construct_semiconductor(T, sc::Dict, input_units::NamedTuple, transformations = missing)
    Semiconductor{T}(sc, input_units, transformations)
end

function construct_passive(T, pass::Dict, input_units::NamedTuple, transformations = missing)
    Passive{T}(pass, input_units, transformations)
end

function construct_contact(T, contact::Dict, input_units::NamedTuple, transformations = missing)
    Contact{T}(contact, input_units, transformations)
end

function construct_virtual_volume(T, pass::Dict, input_units::NamedTuple, transformations = missing)
    construct_virtual_volume(T, pass, input_units, Val{Symbol(pass["model"])}, transformations)
end
function construct_virtual_volume(T, pass::Dict, input_units::NamedTuple, ::Type{Val{:dead}}, transformations = missing)
    DeadVolume{T}(pass, input_units, transformations)
end
function construct_virtual_volume(T, pass::Dict, input_units::NamedTuple, ::Type{Val{:arbitrary}}, transformations = missing)
    ArbitraryDriftModificationVolume{T}(pass, input_units, transformations)
end


function get_world_limits_from_objects(::Type{Cylindrical}, s::Semiconductor{T}, c::Vector{Contact{T}}, p::Vector{Passive{T}}) where {T <: SSDFloat}
    ax1l::T, ax1r::T, ax2l::T, ax2r::T, ax3l::T, ax3r::T = 0, 1, 0, 1, 0, 1
    imps_1::Vector{T} = []
    imps_3::Vector{T} = []
    for objects in [c, p]
        for object in objects
            for posgeo in object.geometry_positive
                append!(imps_1, get_important_points( posgeo, Val{:r}()))
                append!(imps_3, get_important_points( posgeo, Val{:z}()))
            end
        end
    end
    unique!(sort!(imps_1))
    unique!(sort!(imps_3))
    if length(imps_1) > 1
        ax1l = minimum(imps_1)
        ax1r = maximum(imps_1)
    elseif length(imps_1) == 1
        ax1l = minimum(imps_1)
        ax1r = maximum(imps_1) + 1
    end
    if length(imps_3) > 1
        ax3l = minimum(imps_3)
        ax3r = maximum(imps_3)
    elseif length(imps_3) == 1
        ax3l = minimum(imps_3)
        ax3r = maximum(imps_3) + 1
    end
    return ax1l, ax1r, ax2l, ax2r, ax3l, ax3r
end
function get_world_limits_from_objects(::Type{Cartesian}, s::Semiconductor{T}, c::Vector{Contact{T}}, p::Vector{Passive{T}}) where {T <: SSDFloat}
    ax1l::T, ax1r::T, ax2l::T, ax2r::T, ax3l::T, ax3r::T = 0, 1, 0, 1, 0, 1
    imps_1::Vector{T} = []
    imps_2::Vector{T} = []
    imps_3::Vector{T} = []
    for objects in [c, p]
        for object in objects
            for posgeo in object.geometry_positive
                append!(imps_1, get_important_points( posgeo, Val{:x}()))
                append!(imps_2, get_important_points( posgeo, Val{:y}()))
                append!(imps_3, get_important_points( posgeo, Val{:z}()))
            end
        end
    end
    unique!(sort!(imps_1))
    unique!(sort!(imps_2))
    unique!(sort!(imps_3))
    if length(imps_1) > 1
        ax1l = minimum(imps_1)
        ax1r = maximum(imps_1)
    elseif length(imps_1) == 1
        ax1l = minimum(imps_1)
        ax1r = maximum(imps_1) + 1
    end
    if length(imps_2) > 1
        ax2l = minimum(imps_2)
        ax2r = maximum(imps_2)
    elseif length(imps_2) == 1
        ax2l = minimum(imps_2)
        ax2r = maximum(imps_2) + 1
    end
    if length(imps_3) > 1
        ax3l = minimum(imps_3)
        ax3r = maximum(imps_3)
    elseif length(imps_3) == 1
        ax3l = minimum(imps_3)
        ax3r = maximum(imps_3) + 1
    end
    return ax1l, ax1r, ax2l, ax2r, ax3l, ax3r
end

function SolidStateDetector{T}(config_file::Dict, input_units::NamedTuple) where {T <: SSDFloat}
    if haskey(config_file, "detectors")
        config_detector = config_file["detectors"][1] # still only one detector
        
        transformation_keys = filter(k -> k in ("translate", "rotate"), keys(config_detector))
        transformations = broadcast(t -> parse_CSG_transformation(T, config_detector, CSG_dict[t], input_units), transformation_keys)
        if isempty(transformations) transformations = missing end

        @assert haskey(config_detector, "bulk") "Each detector needs an entry `bulk`. Please define the bulk."     
        semiconductor = construct_semiconductor(T, config_detector["bulk"], input_units, transformations)

        @assert haskey(config_detector, "contacts") "Each detector needs at least two contacts. Please define the them in the configuration file."                    
        contacts = broadcast(c -> construct_contact(T, c, input_units, transformations), config_detector["contacts"])
        
        virtual_drift_volumes = if haskey(config_detector, "virtual_drift_volumes")  
            broadcast(v -> construct_virtual_volume(T, v, input_units, transformations), config_detector["virtual_drift_volumes"]) 
        else
            missing
        end
    end
    passives = if haskey(config_file, "surroundings")
        broadcast(p -> construct_passive(T, p, input_units, transformations), config_file["surroundings"])
    else
        missing
    end

    name = haskey(config_file, "name") ? config_file["name"] : "NoNameDetector"
    SolidStateDetector{T}( name, semiconductor, contacts, passives, virtual_drift_volumes )
end

function SolidStateDetector(parsed_dict::Dict)
    SolidStateDetector{Float32}(parsed_dict)
end

function Base.sort!(v::AbstractVector{<:AbstractGeometry})
    hierarchies::Vector{Int} = map(x->x.hierarchy,v)
    v_result::typeof(v) = []
    for idx in unique!(sort!(hierarchies))
        push!(v_result,filter(x->x.hierarchy == hierarchies[idx],v)...)
    end
    return v_result
end

function in(pt::AbstractCoordinatePoint{T}, c::SolidStateDetector{T})::Bool where T
    in(pt,c.semiconductor) || reduce((x,contact) -> x || in(pt,contact), c.contacts, init = false)
end

function println(io::IO, d::SolidStateDetector{T}) where {T <: SSDFloat}
    println("________"*d.name*"________\n")
    println("---General Properties---")
    println("- Precision type: $(T)")
    println()
    println("\t_____Semiconductor_____\n")
    println(d.semiconductor)
    println()
    println("# Contacts: $(length(d.contacts))")
    if length(d.contacts)<=5
        for c in d.contacts
            println(c)
        end
    end
    println()
    println("# Passives: $(length(d.passives))")
    if length(d.passives)<=5
        for p in d.passives
            # println(c)
        end
    end
end

function show(io::IO, d::SolidStateDetector{T}) where {T <: SSDFloat} println(io, d) end
function print(io::IO, d::SolidStateDetector{T}) where {T <: SSDFloat} println(io, d) end
function show(io::IO,::MIME"text/plain", d::SolidStateDetector) where {T <: SSDFloat}
    show(io, d)
end


# ToDo: Test it
function generate_random_startpositions(d::SolidStateDetector{T}, n::Int, Volume::NamedTuple=bounding_box(d), rng::AbstractRNG = MersenneTwister(), min_dist_from_boundary = 0.0001) where T
    delta = T(min_dist_from_boundary)
    n_filled::Int = 0
    positions = Vector{CartesianPoint{T}}(undef,n)
    while n_filled < n
        sample=CylindricalPoint{T}(rand(rng,Volume[:r_range].left:0.00001:Volume[:r_range].right),rand(rng,Volume[:φ_range].left:0.00001:Volume[:φ_range].right),rand(rng,Volume[:z_range].left:0.00001:Volume[:z_range].right))
        if !(sample in d.contacts) && in(sample,d) && in(CylindricalPoint{T}(sample.r+delta,sample.φ,sample.z),d) && in(CylindricalPoint{T}(sample.r-delta,sample.φ,sample.z),d) && in(CylindricalPoint{T}(sample.r,sample.φ,sample.z+delta),d) && in(CylindricalPoint{T}(sample.r,sample.φ,sample.z-delta),d)
            n_filled += 1
            positions[n_filled]=CartesianPoint(sample)
        end
    end
    positions
end

function sample(obj::AbstractObject{T}, g::Grid{T,3,Cartesian,AT})::Vector{CartesianPoint{T}} where {T,AT}
    g_mid::CartesianTicksTuple{T} = ( x = _get_mid_ticks(g[1].ticks), y = _get_mid_ticks(g[2].ticks), z = _get_mid_ticks(g[3].ticks))
    #filter!(p -> p in obj, 
    samples::Vector{CartesianPoint{T}} = vcat([sample(surf, g_mid) for surf in obj.decomposed_surfaces]...)
    unique!(samples)
    #)
end

function sample(obj::AbstractObject{T}, g::Grid{T,3,Cylindrical,AT})::Vector{CylindricalPoint{T}} where {T,AT}
    g_mid::CylindricalTicksTuple{T} = ( r = _get_mid_ticks(g[1].ticks), φ = _get_mid_ticks(g[2].ticks), z = _get_mid_ticks(g[3].ticks))
    #filter!(p -> p in obj, 
    samples::Vector{CylindricalPoint{T}} = vcat([sample(surf, g_mid) for surf in obj.decomposed_surfaces]...)
    unique!(samples)
    #)
end

# TODO: if mid ticks are needed: uncomment. If not: remove
@inline function _get_mid_ticks(gr::Vector{T})::Vector{T} where {T}
    gr
    # gr_mid::Vector{T} = Vector{T}(undef, 2 * length(gr) - 1)
    # gr_mid[1:2:end] = gr
    # gr_mid[2:2:end] = midpoints(gr)
    # gr_mid
end

function paint_object(object::AbstractObject, grid::Grid{T, 3, S})::Vector{NTuple{3,Int}} where {T <: SSDFloat, S}
    samples = sample(object, grid)
    object_gridpoints::Vector{NTuple{3,Int}} = [find_closest_gridpoint(sample_point, grid) for sample_point in samples]
    unique!(object_gridpoints)
end
