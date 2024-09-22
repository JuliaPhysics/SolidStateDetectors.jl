"""
    struct SolidStateDetector{T,SC,CT,PT,VDM} <: AbstractConfig{T}

Struct to describe all parts of a solid state detector, i.e.
the [`Semiconductor`](@ref), a set of [`Contact`](@ref) and (optionally) [`Passive`](@ref) and virtual drift volumes.

The properties of the parts (charge densities, fixed potentials, relative permittivity
of the materials) will be used as input to the calculation of [`ElectricPotential`](@ref) and 
[`WeightingPotential`](@ref) in the [`Simulation`](@ref).

## Parametric types
* `T`: Precision type.
* `SC`: Type of the `semiconductor`.
* `CT`: Type of the `contacts`.
* `PT`: Type of the `passives`.
* `VDM`: Type of the `virtual_drift_volumes`.

## Fields
* `name::String`: Name of the detector.
* `semiconductor::SC`: [`Semiconductor`](@ref) of the detector. 
* `contacts::CT`: Vector of [`Contact`](@ref) of the detector. 
* `passives::PT`: Vector of [`Passive`](@ref) objects, e.g. holding structures around the detector. 
* `virtual_drift_volumes::VDM`: Vector of virtual drift volumes in which the drift can be modulated
    by user-defined methods for `modulate_driftvector`, e.g. [`DeadVolume`](@ref).

See also [`Semiconductor`](@ref), [`Contact`](@ref), [`Passive`](@ref) and [`DeadVolume`](@ref).
"""
struct SolidStateDetector{T,SC,CT,PT,VDM} <: AbstractConfig{T}
    name::AbstractString # optional
    semiconductor::SC
    contacts::CT
    passives::PT
    virtual_drift_volumes::VDM
    
    SolidStateDetector{T}(n::AbstractString,s::SC,c::C,p::P,v::VDM) where {T,SC,C,P,VDM}= new{T,SC,C,P,VDM}(n,s,c,p,v)
end

function SolidStateDetector(det::SolidStateDetector{T}, impurity_density::AbstractImpurityDensity{T}) where {T <: SSDFloat}
    sc = Semiconductor(det.semiconductor, impurity_density)
    SolidStateDetector{T}(
        det.name, sc, det.contacts, det.passives, det.virtual_drift_volumes    
    )
end
function SolidStateDetector(det::SolidStateDetector{T}, charge_drift_model::AbstractChargeDriftModel{T}) where {T <: SSDFloat}
    sc = Semiconductor(det.semiconductor, charge_drift_model)
    SolidStateDetector{T}(
        det.name, sc, det.contacts, det.passives, det.virtual_drift_volumes    
    )
end
function SolidStateDetector(det::SolidStateDetector{T}, charge_trapping_model::AbstractChargeTrappingModel{T}) where {T <: SSDFloat}
    sc = Semiconductor(det.semiconductor, charge_trapping_model)
    SolidStateDetector{T}(
        det.name, sc, det.contacts, det.passives, det.virtual_drift_volumes    
    )
end
function SolidStateDetector(det::SolidStateDetector{T}; contact_id::Int, contact_potential::Real) where {T <: SSDFloat}
    oc = det.contacts[contact_id]
    nc = Contact(T(contact_potential), oc.material, oc.id, oc.name, oc.geometry )
    contacts = [c.id == contact_id ? nc : c for c in det.contacts]
    SolidStateDetector{T}( det.name, det.semiconductor, contacts, det.passives, det.virtual_drift_volumes )
end

get_precision_type(::SolidStateDetector{T}) where {T} = T


function get_world_limits_from_objects(::Type{Cylindrical}, det::SolidStateDetector{T}) where {T <: SSDFloat}
    ax1l::T, ax1r::T, ax2l::T, ax2r::T, ax3l::T, ax3r::T = 0, 1, 0, 1, 0, 1
    t::Array{T, 2} = hcat(vcat(Vector{CylindricalPoint{T}}.(ConstructiveSolidGeometry.extreme_points.(ConstructiveSolidGeometry.surfaces(det.semiconductor.geometry)))...)...)
    (ax1l, ax1r), (ax3l, ax3r) = broadcast(i -> extrema(t[i,:]), 1:2:3)
    for c in det.contacts
        t = hcat([ax1l ax1r; ax2l ax2r; ax3l ax3r], hcat.(Vector{CylindricalPoint{T}}.(ConstructiveSolidGeometry.extreme_points.(ConstructiveSolidGeometry.surfaces(c.geometry)))...)...)
        (ax1l, ax1r), (ax3l, ax3r) = broadcast(i -> extrema(t[i,:]), 1:2:3)
    end
    if !ismissing(det.passives)
        for p in det.passives
            t = hcat([ax1l ax1r; ax2l ax2r; ax3l ax3r], hcat.(Vector{CylindricalPoint{T}}.(ConstructiveSolidGeometry.extreme_points.(ConstructiveSolidGeometry.surfaces(p.geometry)))...)...)
            (ax1l, ax1r), (ax3l, ax3r) = broadcast(i -> extrema(t[i,:]), 1:2:3)
        end
    end
    return ax1l, ax1r, ax2l, ax2r, ax3l, ax3r
end

function get_world_limits_from_objects(::Type{Cartesian}, det::SolidStateDetector{T}) where {T <: SSDFloat}
    ax1l::T, ax1r::T, ax2l::T, ax2r::T, ax3l::T, ax3r::T = 0, 1, 0, 1, 0, 1
    t::Array{T, 2} = hcat(vcat(ConstructiveSolidGeometry.extreme_points.(ConstructiveSolidGeometry.surfaces(det.semiconductor.geometry))...)...)
    (ax1l, ax1r), (ax2l, ax2r), (ax3l, ax3r) = broadcast(i -> extrema(t[i,:]), 1:3)
    for c in det.contacts
        t = hcat([ax1l ax1r; ax2l ax2r; ax3l ax3r], hcat.(ConstructiveSolidGeometry.extreme_points.(ConstructiveSolidGeometry.surfaces(c.geometry))...)...)
        (ax1l, ax1r), (ax2l, ax2r), (ax3l, ax3r) = broadcast(i -> extrema(t[i,:]), 1:3)
    end
    if !ismissing(det.passives)
        for p in det.passives
            t = hcat([ax1l ax1r; ax2l ax2r; ax3l ax3r], hcat.(ConstructiveSolidGeometry.extreme_points.(ConstructiveSolidGeometry.surfaces(p.geometry))...)...)
            (ax1l, ax1r), (ax2l, ax2r), (ax3l, ax3r) = broadcast(i -> extrema(t[i,:]), 1:3)
        end
    end
    return ax1l, ax1r, ax2l, ax2r, ax3l, ax3r
end

function SolidStateDetector{T}(config_file::AbstractDict, input_units::NamedTuple) where {T <: SSDFloat}
    
    @assert !haskey(config_file, "objects") "Configuration file deprecation.\n
        The configuration file format was updated in v0.6.0.
        However, this configuration file still seems to be in the old format.\n
        To update your configuration file to the new format (v0.6.0 and newer),
        open a new Julia session and load the following file:\n
        \tinclude(\"<path_to_SolidStateDetectors.jl>/test/update_config_files.jl\")\n
        Afterwards, run\n
        \tupdate_config_file(\"<path_to_configuration_file>\")\n
        This method returns the file name of the updated configuration file.
        Please close the Julia session after updating the configuration files, as some
        parsing methods are overridden with old methods.\n
        Note that if your old geometries were defined using `difference`, you might
        need to increase the dimensions of the subtracted geometry as it is now treated as
        open primitive. Please have a look at the documentation:\n
        \thttps://juliaphysics.github.io/SolidStateDetectors.jl/stable/\n"
        
    @assert haskey(config_file, "detectors") "Config file needs an entry `detectors` that defines the detector(s)."
    config_detector = config_file["detectors"][1] # still only one detector

    transformations = parse_CSG_transformation(T, config_detector, input_units)
    
    @assert haskey(config_detector, "semiconductor") "Each detector needs an entry `semiconductor`. Please define the semiconductor."     
    semiconductor = Semiconductor{T}(config_detector["semiconductor"], input_units, transformations)

    @assert haskey(config_detector, "contacts") "Each detector needs at least two contacts. Please define them in the configuration file."                    
    contacts = broadcast(c -> Contact{T}(c, input_units, transformations), config_detector["contacts"])
    
    # SolidStateDetectors.jl does not allow for arbitrary contact IDs yet (issue #288)
    # They need to be in order (1, 2, ... , N), so throw an error if this is not the case.
    if !all(getfield.(contacts, :id) .== eachindex(contacts))
        ArgumentError("SolidStateDetectors.jl only supports contact IDs that are in order.\n
            Please set the ID of the first contact to 1, the ID of the second contact to 2, etc.")
    end
    
    passives = Passive{T}[]
    if haskey(config_detector, "passives") # "passives" as entry of "detectors"
        append!(passives, broadcast(p -> Passive{T}(p, input_units, transformations), config_detector["passives"]))
    end
    if haskey(config_file, "surroundings") # "surroundings" as entry in the configuration file
        append!(passives, broadcast(p -> Passive{T}(p, input_units, parse_CSG_transformation(T, p, input_units)), config_file["surroundings"]))
    end
    if isempty(passives) 
        passives = missing 
    end
        
    virtual_drift_volumes = if haskey(config_detector, "virtual_drift_volumes")  
        broadcast(v -> construct_virtual_volume(T, v, input_units, transformations), config_detector["virtual_drift_volumes"]) 
    else
        missing
    end

    name = haskey(config_file, "name") ? config_file["name"] : "NoNameDetector"
    SolidStateDetector{T}( name, semiconductor, contacts, passives, virtual_drift_volumes )
end

function SolidStateDetector(dict::AbstractDict)
    SolidStateDetector{Float32}(dict)
end

function in(pt::AbstractCoordinatePoint{T}, det::SolidStateDetector{T})::Bool where T
    in(pt, det.semiconductor) || reduce((x, contact) -> x || in(pt, contact), det.contacts, init = false)
end

function println(io::IO, det::SolidStateDetector{T}) where {T <: SSDFloat}
    println("________"*det.name*"________\n")
    println("---General Properties---")
    println("- Precision type: $(T)")
    println()
    println("\t_____Semiconductor_____\n")
    println(det.semiconductor)
    println()
    println("# Contacts: $(length(det.contacts))")
    if length(det.contacts)<=5
        for contact in det.contacts
            println(contact)
        end
    end
    if !ismissing(det.passives)
        println()
        println("# Passives: $(length(det.passives))")
        if length(det.passives) <= 5
            for passive in det.passives
                println(passive)
            end
        end
    end
end

function show(io::IO, det::SolidStateDetector{T}) where {T <: SSDFloat} println(io, det) end
function print(io::IO, det::SolidStateDetector{T}) where {T <: SSDFloat} println(io, det) end
function show(io::IO,::MIME"text/plain", det::SolidStateDetector{T}) where {T <: SSDFloat} show(io, det) end

function determine_bias_voltage_contact_id(det::SolidStateDetector{T}) where {T <: SSDFloat}
    contact_potentials = T[c.potential for c in det.contacts]
    inds = findall(!iszero, contact_potentials)
    @assert length(inds) == 1 "Could not determine contact at which the bias voltage is applied as multiple contacts have non-zero contact potentials."
    inds[1]
end
