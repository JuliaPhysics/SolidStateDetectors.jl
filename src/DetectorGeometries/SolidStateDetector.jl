"""
    mutable struct SolidStateDetector{T <: SSDFloat, CS} <: AbstractConfig{T}

CS: Coordinate System: -> :Cartesian / :Cylindrical
"""
mutable struct SolidStateDetector{T <: SSDFloat, CS} <: AbstractConfig{T}
    name::String  # optional
    class::Symbol # optional

    material_detector::NamedTuple
    material_environment::NamedTuple
    cyclic::T # optional
    mirror_symmetry_φ::Bool # optional

    geometry_unit::Unitful.Units
    geometry_unit_factor::Real # optional, geometry_unit is enough
    contacts_geometry_unit::Unitful.Units # optional, why not geometry_unit

    world::AbstractGeometry{T}

    crystal_geometry::AbstractGeometry{T}

    bulk_type::Symbol
    charge_density_model::AbstractChargeDensityModel{T}

    contacts::Vector{AbstractContact{T}}

    external_parts::Vector{AbstractContact{T}}
    geometry_external_positive::Vector{AbstractGeometry{T}} # ?
    geometry_external_negative::Vector{AbstractGeometry{T}} # ?

    SolidStateDetector{T, CS}() where {T <: SSDFloat, CS} = new{T, CS}()
end

get_precision_type(d::SolidStateDetector{T}) where {T} = T
get_coordinate_system(d::SolidStateDetector{T, CS}) where {T, CS} = CS

function SolidStateDetector{T}(config_file::Dict)::SolidStateDetector{T} where{T <: SSDFloat}
    c = if Symbol(config_file["class"]) == :CGD
        SolidStateDetector{T, :Cartesian}()
    else
        SolidStateDetector{T, :Cylindrical}()
    end
    c.class = Symbol(config_file["class"])
    c.name = config_file["name"]
    if c.class != :CGD
        c.cyclic = T(deg2rad(config_file["cyclic"]))
        c.mirror_symmetry_φ = config_file["mirror_symmetry_phi"] == "true"
    end

    c.material_environment = material_properties[materials[config_file["geometry"]["world"]["material"]]]
    c.material_detector = material_properties[materials[config_file["geometry"]["crystal"]["material"]]]

    c.geometry_unit = unit_conversion[config_file["geometry"]["unit"]]
    c.world = Geometry(T, config_file["geometry"]["world"]["geometry"], c.geometry_unit)[2][1]

    haskey(config_file["geometry"],"external") ? c.external_parts = Contact{T,:E}[ Contact{T,:E}( ep_dict, c.geometry_unit) for ep_dict in config_file["geometry"]["external"]] : c.external_parts = []

    c.crystal_geometry , geometry_positive, geometry_negative = Geometry(T, config_file["geometry"]["crystal"]["geometry"], c.geometry_unit)

    c.bulk_type = bulk_types[config_file["bulk_type"]]

    c.charge_density_model = ChargeDensityModel(T, config_file["charge_density_model"])

    c.contacts_geometry_unit = unit_conversion[config_file["contacts"]["unit"]]
    haskey(config_file["contacts"], "p") ? p_contacts = Contact{T, :P}[ Contact{T, :P}( contact_dict, c.geometry_unit ) for contact_dict in config_file["contacts"]["p"] ] : nothing
    haskey(config_file["contacts"], "n") ? n_contacts = Contact{T, :N}[ Contact{T, :N}( contact_dict, c.geometry_unit ) for contact_dict in config_file["contacts"]["n"] ] : nothing
    c.contacts = vcat(p_contacts, n_contacts)

    return c
end

function SolidStateDetector(parsed_dict::Dict)
    SolidStateDetector{Float32}(parsed_dict)
end

function Base.sort!(v::AbstractVector{<:AbstractGeometry})
    hierarchies::Vector{Int} = map(x->x.hierarchy,v)
    v_result::typeof(v) = []
    for idx in sort!(unique!(hierarchies))
        push!(v_result,filter(x->x.hierarchy == hierarchies[idx],v)...)
    end
    return v_result
end

function SolidStateDetector{T}(parsed_dict::Dict) where T
    SolidStateDetector{T}(parsed_dict)
end

function contains(c::SolidStateDetector, point::AbstractCoordinatePoint{T,3})::Bool where T
    for contact in c.contacts
        if point in contact
            return true
        end
    end
    if point in c.crystal_geometry
        return true
    else
        return false
    end
end

function println(io::IO, d::SolidStateDetector{T}) where {T <: SSDFloat}
    println("________"*d.name*"________\n")
    println("Class: ",d.class)
    println("---General Properties---")
    println("Detector Material: \t $(d.material_detector.name)")
    println("Environment Material: \t $(d.material_environment.name)")
    println("Bulk type: \t\t $(d.bulk_type)")
    # println("Core Bias Voltage: \t $(d.segment_bias_voltages[1]) V")
    # println("Mantle Bias Voltage: \t $(d.segment_bias_voltages[2]) V\n")
end

function show(io::IO, d::SolidStateDetector{T}) where {T <: SSDFloat} println(d) end
function print(io::IO, d::SolidStateDetector{T}) where {T <: SSDFloat} println(d) end
function display(io::IO, d::SolidStateDetector{T} ) where {T <: SSDFloat} println(d) end
function show(io::IO,::MIME"text/plain", d::SolidStateDetector) where {T <: SSDFloat}
    show(io, d)
end
