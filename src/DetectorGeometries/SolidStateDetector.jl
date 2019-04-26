"""
    mutable struct SolidStateDetector{T <: SSDFloat, CS} <: AbstractConfig{T}

CS: Coordinate System: -> :Cartesian / :Cylindrical
"""
mutable struct SolidStateDetector{T <: SSDFloat, CS} <: AbstractConfig{T}
    name::String  # optional
    inputunits::Dict{String,Unitful.Units}
    world::AbstractVolumePrimitive
    grid_type::Symbol
    mirror_symmetry_φ::Bool # optional
    cyclic::T
    medium::NamedTuple

    semiconductors::Vector{Semiconductor{T}}
    contacts::Vector{Contact{T}}
    passives::Vector{Passive{T}}

    SolidStateDetector{T, CS}() where {T <: SSDFloat, CS} = new{T, CS}()
end

get_precision_type(d::SolidStateDetector{T}) where {T} = T
get_coordinate_system(d::SolidStateDetector{T, CS}) where {T, CS} = CS

function construct_units(dict::Dict)::Dict{String,Unitful.Units}
    result_dict::Dict{String,Unitful.Units} = Dict()
    haskey(dict,"length") ? result_dict["length"] = unit_conversion[dict["length"]] : result_dict["length"] = u"mm"
    haskey(dict,"angle") ? result_dict["angle"] = unit_conversion[dict["angle"]] : result_dict["angle"] = u"rad"
    haskey(dict,"potential") ? result_dict["potential"] = unit_conversion[dict["potential"]] : result_dict["potential"] = u"V"
    haskey(dict,"temperature") ? result_dict["temperature"] = unit_conversion[dict["temperature"]] : result_dict["temperature"] = u"K"
    result_dict
end

function construct_world(T, dict::Dict, inputunit_dict::Dict{String, Unitful.Units})::Tuple{Symbol,AbstractVolumePrimitive}
    if dict["coordinates"] == "Cylindrical"
        vol = Tube{T}(
        Interval(geom_round(ustrip(uconvert(u"m", T(0.0) * inputunit_dict["length"] ))), geom_round(ustrip(uconvert(u"m", T(dict["dimensions"]["r"]) * inputunit_dict["length"])))),
        Interval(geom_round(ustrip(uconvert(u"rad",T(0.0) * inputunit_dict["angle"] ))),geom_round(ustrip(uconvert(u"rad",T(360.0) * inputunit_dict["angle"]))) ),
        Interval(geom_round(ustrip(uconvert(u"m", T(dict["dimensions"]["z"]["from"]) * inputunit_dict["length"] ))), geom_round(ustrip(uconvert(u"m", T(dict["dimensions"]["z"]["to"]) * inputunit_dict["length"])))),
        missing
        )
    elseif dict["coordinates"] == "Cartesian"
        # if typeof(dict["dimensions"]["x"]) <: SSDFloat
        #     xStart, xStop = geom_round(translate[1] + T(0.0)), geom_round(translate[1] + ustrip(uconvert(u"m", T(dict["dimensions"]["x"]) * inputunit_dict["length"])))
        # else
        #     xStart, xStop = geom_round(translate[1] + ustrip(uconvert(u"m", T(dict["dimensions"]["x"]["from"]) * inputunit_dict["length"] ))), geom_round(translate[1] + ustrip(uconvert(u"m", T(dict["dimensions"]["x"]["to"]) * inputunit_dict["length"])))
        # end
        # if typeof(dict["dimensions"]["y"]) <: SSDFloat
        #     yStart, yStop = geom_round(translate[1] + T(0.0)), geom_round(translate[1] + ustrip(uconvert(u"m", T(dict["dimensions"]["y"]) * inputunit_dict["length"])))
        # else
        #     yStart, yStop = geom_round(translate[1] + ustrip(uconvert(u"m", T(dict["dimensions"]["y"]["from"]) * inputunit_dict["length"] ))), geom_round(translate[1] + ustrip(uconvert(u"m", T(dict["dimensions"]["y"]["to"]) * inputunit_dict["length"])))
        # end
        # if typeof(dict["dimensions"]["z"]) <: SSDFloat
        #     zStart, zStop = geom_round(translate[1] + T(0.0)), geom_round(translate[1] + ustrip(uconvert(u"m", T(dict["dimensions"]["z"]) * inputunit_dict["length"])))
        # else
        #     zStart, zStop = geom_round(translate[1] + ustrip(uconvert(u"m", T(dict["dimensions"]["z"]["from"]) * inputunit_dict["length"] ))), geom_round(translate[1] + ustrip(uconvert(u"m", T(dict["dimensions"]["z"]["to"]) * inputunit_dict["length"])))
        # end
        # translate = missing
        # vol= CartesianBox3D{T}(
        #     (xStart, xStop),
        #     (yStart, yStop),
        #     (zStart, zStop),
        #     translate )
        vol = CartesianBox3D{T}(dict["dimensions"], inputunit_dict)
    else
        @warn "Gridtype must be 'Cylindrical' or 'Cartesian'"
    end
    return Symbol(dict["coordinates"]), vol
end

function construct_semiconductor(T, sc::Dict, inputunit_dict::Dict{String, Unitful.Units})
    Semiconductor{T}(sc, inputunit_dict)
end

function construct_passive(T, pass::Dict, inputunit_dict::Dict{String, Unitful.Units})
    Passive{T}(pass, inputunit_dict)
end

function construct_contact(T, contact::Dict, inputunit_dict::Dict{String, Unitful.Units})
    Contact{T}(contact, inputunit_dict)
end

function construct_objects(T, objects::Vector, semiconductors, contacts, passives, inputunit_dict)::Nothing
    for obj in objects
        if obj["class"] == "Semiconductor"
            push!(semiconductors, construct_semiconductor(T, obj, inputunit_dict))
        elseif obj["class"] == "Contact"
            push!(contacts, construct_contact(T, obj, inputunit_dict))
        elseif obj["class"] == "Passive"
            push!(passives, construct_passive(T, obj, inputunit_dict))
        else
            @warn "please spcify the calss to bei either a 'Semiconductor', a 'Contact', or 'Passive'"
        end
    end
    nothing
end

function SolidStateDetector{T}(config_file::Dict)::SolidStateDetector{T} where{T <: SSDFloat}
    grid_type = Symbol(config_file["world"]["grid"]["coordinates"])
    c = SolidStateDetector{T, grid_type}()
    c.name = config_file["name"]
    c.inputunits = construct_units(config_file["world"]["units"])
    c.grid_type, c.world = construct_world(T, config_file["world"]["grid"], c.inputunits)
    c.medium = material_properties[materials[config_file["world"]["medium"]]]
    c.semiconductors, c.contacts, c.passives = [], [], []
    c.cyclic = geom_round(T(ustrip(uconvert(u"rad",config_file["world"]["grid"]["symmetries"]["periodic"]["phi"] * c.inputunits["angle"]))))
    c.mirror_symmetry_φ = false
    construct_objects(T, config_file["world"]["objects"], c.semiconductors, c.contacts, c.passives, c.inputunits)

    # c = if Symbol(config_file["class"]) == :CGD
    #     SolidStateDetector{T, :Cartesian}()
    # else
    #     SolidStateDetector{T, :Cylindrical}()
    # end
    # c.class = Symbol(config_file["class"])
    # c.name = config_file["name"]
    # if c.class != :CGD
    #     c.cyclic = T(deg2rad(config_file["cyclic"]))
    #     c.mirror_symmetry_φ = config_file["mirror_symmetry_phi"] == "true"
    # end
    #
    # c.material_environment = material_properties[materials[config_file["geometry"]["world"]["material"]]]
    # c.material_detector = material_properties[materials[config_file["geometry"]["crystal"]["material"]]]
    #
    # c.geometry_unit = unit_conversion[config_file["geometry"]["unit"]]
    # c.world = Geometry(T, config_file["geometry"]["world"]["geometry"], c.geometry_unit)[2][1]
    #
    # haskey(config_file["geometry"],"external") ? c.external_parts = Contact{T,:E}[ Contact{T,:E}( ep_dict, c.geometry_unit) for ep_dict in config_file["geometry"]["external"]] : c.external_parts = []
    #
    # c.crystal_geometry , geometry_positive, geometry_negative = Geometry(T, config_file["geometry"]["crystal"]["geometry"], c.geometry_unit)
    #
    # c.bulk_type = bulk_types[config_file["bulk_type"]]
    #
    # c.charge_density_model = ChargeDensityModel(T, config_file["charge_density_model"])
    #
    # c.contacts_geometry_unit = unit_conversion[config_file["contacts"]["unit"]]
    # haskey(config_file["contacts"], "p") ? p_contacts = Contact{T, :P}[ Contact{T, :P}( contact_dict, c.geometry_unit ) for contact_dict in config_file["contacts"]["p"] ] : nothing
    # haskey(config_file["contacts"], "n") ? n_contacts = Contact{T, :N}[ Contact{T, :N}( contact_dict, c.geometry_unit ) for contact_dict in config_file["contacts"]["n"] ] : nothing
    # c.contacts = vcat(p_contacts, n_contacts)

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
    for sc in c.semiconductors
        if point in sc
            return true
        end
    end
    return false
end

function println(io::IO, d::SolidStateDetector{T}) where {T <: SSDFloat}
    println("________"*d.name*"________\n")
    # println("Class: ",d.class)
    println("---General Properties---")
    # for d in semiconductors
    #     println("Detector Material: \t $(d.material_detector.name)")
    #     println("Environment Material: \t $(d.material_environment.name)")
    #     println("Bulk type: \t\t $(d.bulk_type)")
    # end

    # println("Core Bias Voltage: \t $(d.segment_bias_voltages[1]) V")
    # println("Mantle Bias Voltage: \t $(d.segment_bias_voltages[2]) V\n")
end

function show(io::IO, d::SolidStateDetector{T}) where {T <: SSDFloat} println(d) end
function print(io::IO, d::SolidStateDetector{T}) where {T <: SSDFloat} println(d) end
function display(io::IO, d::SolidStateDetector{T} ) where {T <: SSDFloat} println(d) end
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
        if !(sample in d.contacts) && contains(d,sample) && contains(d,CylindricalPoint{T}(sample.r+delta,sample.φ,sample.z))&& contains(d,CylindricalPoint{T}(sample.r-delta,sample.φ,sample.z))&& contains(d,CylindricalPoint{T}(sample.r,sample.φ,sample.z+delta))&& contains(d,CylindricalPoint{T}(sample.r,sample.φ,sample.z-delta))
            n_filled += 1
            positions[n_filled]=CartesianPoint(sample)
        end
    end
    positions
end
