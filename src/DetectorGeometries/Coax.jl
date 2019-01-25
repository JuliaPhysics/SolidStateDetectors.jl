mutable struct Coax{T<:AbstractFloat} <: SolidStateDetector{T}
    name::String
    material_detector::NamedTuple
    material_environment::NamedTuple
    cyclic::T
    mirror_symmetry_θ::Bool
    path_to_configfile::String
    #### Importet Values from JSON file
    ### Geometry
    geometry_unit_factor::T
    bulk_type::Symbol
    crystal_length::T
    crystal_radius::T

    borehole_length::T
    borehole_radius::T
    borehole_modulation::Bool
    borehole_ModulationFunction::Function
    borehole_segment_idx::Int
    borehole_top_segment_idx::Int
    borehole_bot_segment_idx::Int

    groove_endplate::String
    groove_depth::T
    groove_rInner::T
    groove_width::T

    taper_inner_bot_length::T
    taper_inner_bot_rOuter::T
    taper_inner_top_length::T
    taper_inner_top_rOuter::T

    taper_outer_bot_length::T
    taper_outer_bot_rInner::T
    taper_outer_top_length::T
    taper_outer_top_rInner::T

    charge_carrier_density_unit_factor::T
    charge_carrier_density_top::T
    charge_carrier_density_bot::T

    ### Segmentaion
    n_total_contacts::Int
    # n_repetitive_segments::Int
    n_individual_segments::Int

    segmentation_r_ranges::Array{Tuple{T,T},1}
    segmentation_phi_ranges::Array{Tuple{T,T},1}
    segmentation_z_ranges::Array{Tuple{T,T},1}
    segmentation_types::Array{String,1}
    segment_bias_voltages::Array{T,1}

    segmentation_boundaryWidths_horizontal::Array{Tuple{T,T},1}
    segmentation_boundaryWidths_vertical::Array{Tuple{T,T},1}
    segmentation_boundaryMidpoints_radial::Array{T,1}
    segmentation_boundaryMidpoints_vertical::Array{T,1}
    segmentation_boundaryMidpoints_horizontal::Array{T,1}

    #Grouping of Segments to form Channels (if true, as many as there are n_total_contacts)
    custom_grouping::Bool
    grouped_channels::Array{Array{Int,1},1}

    floating_boundary_r_ranges::Array{Tuple{T,T},1}
    floating_boundary_phi_ranges::Array{Tuple{T,T},1}
    floating_boundary_z_ranges::Array{Tuple{T,T},1}
    floating_boundary_types::Array{String,1}
    ##Generator
    # Coax{T}() where {T<:AbstractFloat} = new("Unknown",-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,Shape([-1],[-1]),)
    Coax{T}() where {T<:AbstractFloat} = new{T}()
    # function Coax{T}() where {T<:AbstactFloat}
    # new()
    # end
end

function Coax(mytype::Type{<:AbstractFloat},config_file::Dict)::Coax
    return  Coax{mytype}(config_file::Dict)
end

function Coax{T}(config_file::Dict)::Coax where T<:AbstractFloat
    config_file["class"] != "Coax" ? error() : nothing
    c::Coax = Coax{T}()
    c.name = config_file["name"]
    unit_conversion = Dict{String, T}("nm"=>.1e-9,"um"=>1e-6,"mm"=>1e-3,"cm"=>1e-2,"m"=>1)
    c.geometry_unit_factor  = unit_conversion[config_file["geometry"]["unit"]]
    materials = Dict("HPGe" => :HPGe, "Vacuum" => :Vacuum)
    c.material_detector = material_properties[materials[config_file["materials"]["detector"]]]
    c.material_environment = material_properties[materials[config_file["materials"]["environment"]]]
    bulk_types = Dict{Any, Symbol}( "n" => :ntype,
        "n-type" => :ntype,
        "ntype" => :ntype,
        "p-type" => :ptype,
        "ptype" => :ptype,
        "p" => :ptype  )
    c.bulk_type = bulk_types[ config_file["type"]  ]
    c.cyclic = deg2rad(config_file["cyclic"])
    c.mirror_symmetry_θ = config_file["mirror_symmetry_θ"] == "true"
    c.borehole_modulation=false
    c.crystal_length = round(c.geometry_unit_factor * config_file["geometry"]["crystal"]["length"], sigdigits=5)
    c.crystal_radius = round(c.geometry_unit_factor * config_file["geometry"]["crystal"]["radius"], sigdigits=5)
    c.borehole_length = round(c.geometry_unit_factor * config_file["geometry"]["borehole"]["length"],sigdigits=5)
    c.borehole_radius = round(c.geometry_unit_factor * config_file["geometry"]["borehole"]["radius"],sigdigits=5)
    # c.pc_depth = round(c.geometry_unit_factor * config_file["geometry"]["point_contact"]["depth"],digits=5)
    # c.pc_radius = round(c.geometry_unit_factor * config_file["geometry"]["point_contact"]["radius"],digits=5)
    c.groove_endplate = config_file["geometry"]["groove"]["endplate"]
    c.groove_depth = round(c.geometry_unit_factor * config_file["geometry"]["groove"]["depth"],digits=5)
    c.groove_rInner = round(c.geometry_unit_factor * config_file["geometry"]["groove"]["rInner"],digits=5)
    c.groove_width = round(c.geometry_unit_factor * config_file["geometry"]["groove"]["width"],digits=5)
    c.taper_inner_bot_length = round(c.geometry_unit_factor * config_file["geometry"]["taper"]["inner"]["bot"]["length"],digits=5)
    c.taper_inner_bot_rOuter = round(c.geometry_unit_factor * config_file["geometry"]["taper"]["inner"]["bot"]["rOuter"],digits=5)
    c.taper_inner_top_length = round(c.geometry_unit_factor * config_file["geometry"]["taper"]["inner"]["top"]["length"],digits=5)
    c.taper_inner_top_rOuter = round(c.geometry_unit_factor * config_file["geometry"]["taper"]["inner"]["top"]["rOuter"],digits=5)
    c.taper_outer_bot_length = round(c.geometry_unit_factor * config_file["geometry"]["taper"]["outer"]["bot"]["length"],digits=5)
    c.taper_outer_bot_rInner = round(c.geometry_unit_factor * config_file["geometry"]["taper"]["outer"]["bot"]["rInner"],digits=5)
    c.taper_outer_top_length = round(c.geometry_unit_factor * config_file["geometry"]["taper"]["outer"]["top"]["length"],digits=5)
    c.taper_outer_top_rInner = round(c.geometry_unit_factor * config_file["geometry"]["taper"]["outer"]["top"]["rInner"],digits=5)
    # c.sfc_depth = round(c.geometry_unit_factor * config_file["geometry"]["surface_contact"]["depth"],digits=5)
    c.charge_carrier_density_top = round(config_file["charge_carrier_density"]["top"],digits=5)
    c.charge_carrier_density_bot = round(config_file["charge_carrier_density"]["bot"],digits=5)

    ### Segmentation
    c.n_total_contacts = config_file["segmentation"]["n_contacts_total"]
    # c.n_repetitive_segments = config_file["segmentation"]["n_repetitive_segments"]
    c.n_individual_segments = config_file["segmentation"]["n_individual_segments"]
    c.custom_grouping = config_file["segmentation"]["custom_grouping"]

    # c.segmentation_boundaryWidth_vertical = round(c.geometry_unit_factor *config_file["segmentation"]["repetitive_segment"]["boundaryWidth"]["vertical"],digits=5)
    # c.segmentation_boundaryWidth_horizontal = round(config_file["segmentation"]["repetitive_segment"]["boundaryWidth"]["vertical"]*c.geometry_unit_factor *180/(π*c.crystal_radius), digits=5)
    # construct_segmentation_arrays_from_repetitive_segment(c, config_file)
    c.n_individual_segments > 0 ? construct_segmentation_arrays_for_individual_segments(c,config_file) : nothing
    construct_grouped_channels_array(c, config_file)
    construct_floating_boundary_arrays(c)
    return c
end

function Coax(mytype::Type{<:AbstractFloat},inputfilename::String)::Coax
    dicttext = read(inputfilename, String)
    parsed_json_file = JSON.parse(dicttext)
    return Coax{mytype}(parsed_json_file)
end
function Coax{T}(inputfilename::String)::Coax where T <: AbstractFloat
    dicttext = read(inputfilename, String)
    parsed_json_file = JSON.parse(dicttext)
    return Coax{T}(parsed_json_file)
end

function Coax(inputfilename::String)::Coax
    # parsed_json_file = Dict()
    # open(inputfilename,"r") do f
    #     dicttext::String = readstring(f)
    #     parsed_json_file = JSON.parse(dicttext)
    # end
    dicttext = read(inputfilename, String)
    parsed_json_file = JSON.parse(dicttext)
    return Coax{Float32}(parsed_json_file)
end

function println(c::Coax)
    println("Coaxial Detector:")
    println("++++++++++++++++++++++")
    println("Geometry:")
    println("Crystal outer Measures: Diameter = $(2*c.crystal_radius), Length = $(c.crystal_length)")
    println("Borehole Diameter = $(2*c.borehole_radius)")
    println("Borehole Tapers: top: r = $(c.taper_inner_top_rOuter), z = $(c.crystal_length) to r = $(c.borehole_radius) ,z = $(c.crystal_length-c.taper_inner_top_length)")
    println("Borehole Tapers: top: r = $(c.taper_inner_bot_rOuter), z = 0 to r = $(c.borehole_radius) ,z = $(c.taper_inner_bot_length)")
end
function show(c::Coax)  println(c)   end
function display(c::Coax)  println(c)   end
function print(c::Coax)  println(c)   end


