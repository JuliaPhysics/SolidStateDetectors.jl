"""
    struct SolidStateDetector{T,SC,CT,PT,VDM} <: AbstractConfig{T}

# Fields
- `name::String`: Name of the detector
- `semiconductor::SC`: Semiconductor of the detector. 
- `contacts::CT`: Vector of `Contact`s of the detector. 
- `passives::PT`: Vector of `Passive` objects, e.g. holding structures around the detector. 
- `virtual_drift_volumes::VDM`: Vector of virtual drift volumes in which the drift can be modulated
    by an user-defined method for `modulate_driftvector`.

"""
struct SolidStateDetector{T,SC,CT,PT,VDM} <: AbstractConfig{T}
    name::String  # optional
    semiconductor::SC
    contacts::CT
    passives::PT
    virtual_drift_volumes::VDM
    
    SolidStateDetector{T}(n::AbstractString,s::SC,c::C,p::P,v::VDM) where {T,SC,C,P,VDM}= new{T,SC,C,P,VDM}(n,s,c,p,v)
end

function SolidStateDetector(det::SolidStateDetector{T,SC,CT,PT,VDM}, impurity_density::AbstractImpurityDensity{T}) where {T,SC,CT,PT,VDM}
    sc = Semiconductor(det.semiconductor, impurity_density)
    SolidStateDetector{T}(
        det.name, sc, det.contacts, det.passives, det.virtual_drift_volumes    
    )
end
function SolidStateDetector(det::SolidStateDetector{T,SC,CT,PT,VDM}, chargedriftmodel::AbstractChargeDriftModel{T}) where {T,SC,CT,PT,VDM}
    sc = Semiconductor(det.semiconductor, chargedriftmodel)
    SolidStateDetector{T}(
        det.name, sc, det.contacts, det.passives, det.virtual_drift_volumes    
    )
end
function SolidStateDetector(det::SolidStateDetector{T,SC,CT,PT,VDM}; contact_id::Int, contact_potential::Real) where {T,SC,CT,PT,VDM}
    oc = det.contacts[contact_id]
    nc = Contact(T(contact_potential), oc.material, oc.id, oc.name, oc.geometry )
    contacts = [c.id == contact_id ? nc : c for c in det.contacts]
    SolidStateDetector{T}( det.name, det.semiconductor, contacts, det.passives, det.virtual_drift_volumes )
end

get_precision_type(::SolidStateDetector{T}) where {T} = T

function construct_virtual_volume(T, pass::Dict, input_units::NamedTuple, transformations::Transformations)
    construct_virtual_volume(T, pass, input_units, Val{Symbol(pass["model"])}, transformations)
end
function construct_virtual_volume(T, pass::Dict, input_units::NamedTuple, ::Type{Val{:dead}}, transformations::Transformations)
    DeadVolume{T}(pass, input_units, transformations)
end
function construct_virtual_volume(T, pass::Dict, input_units::NamedTuple, ::Type{Val{:arbitrary}}, transformations::Transformations)
    ArbitraryDriftModificationVolume{T}(pass, input_units, transformations)
end


function get_world_limits_from_objects(::Type{Cylindrical}, ssd::SolidStateDetector{T}) where {T <: SSDFloat}
    ax1l::T, ax1r::T, ax2l::T, ax2r::T, ax3l::T, ax3r::T = 0, 1, 0, 1, 0, 1
    t::Array{T, 2} = hcat(hcat.(Vector{CylindricalPoint{T}}.(ConstructiveSolidGeometry.extreme_points.(ConstructiveSolidGeometry.surfaces(ssd.semiconductor.geometry)))...)...)
    (ax1l, ax1r), (ax3l, ax3r) = broadcast(i -> extrema(t[i,:]), 1:2:3)
    for c in ssd.contacts
        t = hcat([ax1l ax1r; ax2l ax2r; ax3l ax3r], hcat.(Vector{CylindricalPoint{T}}.(ConstructiveSolidGeometry.extreme_points.(ConstructiveSolidGeometry.surfaces(c.geometry)))...)...)
        (ax1l, ax1r), (ax3l, ax3r) = broadcast(i -> extrema(t[i,:]), 1:2:3)
    end
    if !ismissing(ssd.passives)
        for p in ssd.passives
            t = hcat([ax1l ax1r; ax2l ax2r; ax3l ax3r], hcat.(Vector{CylindricalPoint{T}}.(ConstructiveSolidGeometry.extreme_points.(ConstructiveSolidGeometry.surfaces(p.geometry)))...)...)
            (ax1l, ax1r), (ax3l, ax3r) = broadcast(i -> extrema(t[i,:]), 1:2:3)
        end
    end
    return ax1l, ax1r, ax2l, ax2r, ax3l, ax3r
end

function get_world_limits_from_objects(::Type{Cartesian}, ssd::SolidStateDetector{T}) where {T <: SSDFloat}
    ax1l::T, ax1r::T, ax2l::T, ax2r::T, ax3l::T, ax3r::T = 0, 1, 0, 1, 0, 1
    t::Array{T, 2} = hcat(hcat.(ConstructiveSolidGeometry.extreme_points.(ConstructiveSolidGeometry.surfaces(ssd.semiconductor.geometry))...)...)
    (ax1l, ax1r), (ax2l, ax2r), (ax3l, ax3r) = broadcast(i -> extrema(t[i,:]), 1:3)
    for c in ssd.contacts
        t = hcat([ax1l ax1r; ax2l ax2r; ax3l ax3r], hcat.(ConstructiveSolidGeometry.extreme_points.(ConstructiveSolidGeometry.surfaces(c.geometry))...)...)
        (ax1l, ax1r), (ax2l, ax2r), (ax3l, ax3r) = broadcast(i -> extrema(t[i,:]), 1:3)
    end
    if !ismissing(ssd.passives)
        for p in ssd.passives
            t = hcat([ax1l ax1r; ax2l ax2r; ax3l ax3r], hcat.(ConstructiveSolidGeometry.extreme_points.(ConstructiveSolidGeometry.surfaces(p.geometry))...)...)
            (ax1l, ax1r), (ax2l, ax2r), (ax3l, ax3r) = broadcast(i -> extrema(t[i,:]), 1:3)
        end
    end
    return ax1l, ax1r, ax2l, ax2r, ax3l, ax3r
end

function SolidStateDetector{T}(config_file::Dict, input_units::NamedTuple) where {T <: SSDFloat}
    if haskey(config_file, "detectors")
        config_detector = config_file["detectors"][1] # still only one detector

        transformations = parse_CSG_transformation(T, config_detector, input_units)
        
        @assert haskey(config_detector, "semiconductor") "Each detector needs an entry `semiconductor`. Please define the semiconductor."     
        semiconductor = Semiconductor{T}(config_detector["semiconductor"], input_units, transformations)

        @assert haskey(config_detector, "contacts") "Each detector needs at least two contacts. Please define the them in the configuration file."                    
        contacts = broadcast(c -> Contact{T}(c, input_units, transformations), config_detector["contacts"])
        
        virtual_drift_volumes = if haskey(config_detector, "virtual_drift_volumes")  
            broadcast(v -> construct_virtual_volume(T, v, input_units, transformations), config_detector["virtual_drift_volumes"]) 
        else
            missing
        end
    end
    passives = if haskey(config_file, "surroundings")
        config_surroundings = config_file["surroundings"]
        broadcast(p -> Passive{T}(p, input_units, parse_CSG_transformation(T, p, input_units)), config_file["surroundings"])
    else
        missing
    end

    name = haskey(config_file, "name") ? config_file["name"] : "NoNameDetector"
    SolidStateDetector{T}( name, semiconductor, contacts, passives, virtual_drift_volumes )
end

function SolidStateDetector(parsed_dict::Dict)
    SolidStateDetector{Float32}(parsed_dict)
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
    if !ismissing(d.passives)
        println()
        println("# Passives: $(length(d.passives))")
        if length(d.passives) <= 5
            for p in d.passives
                println(p)
            end
        end
    end
end

function show(io::IO, d::SolidStateDetector{T}) where {T <: SSDFloat} println(io, d) end
function print(io::IO, d::SolidStateDetector{T}) where {T <: SSDFloat} println(io, d) end
function show(io::IO,::MIME"text/plain", d::SolidStateDetector) where {T <: SSDFloat} show(io, d) end

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