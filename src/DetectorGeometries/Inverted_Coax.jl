mutable struct InvertedCoax{T<:AbstractFloat} <: SolidStateDetector{T}
    name::String
    material_detector::NamedTuple
    material_environment::NamedTuple
    cyclic::T
    mirror_symmetry_θ::Bool
    #### Importet Values from JSON file
    geometry_unit_factor
    bulk_type::Symbol
    crystal_length::T
    crystal_radius::T
    pc_depth::T
    pc_radius::T
    groove_depth::T
    groove_rInner::T
    groove_width::T
    groove_endplate::String
    borehole_length::T
    borehole_radius::T
    taper_outer_length::T
    taper_outer_rInner::T
    taper_outer_angle::T
    taper_inner_length::T
    taper_inner_rOuter::T
    taper_inner_angle::T
    borehole_modulation::Bool
    borehole_ModulationFunction::Function
    borehole_segment_idx::Int
    borehole_top_segment_idx::Int
    borehole_bot_segment_idx::Int

    volumes::Array{Volume,1}

    charge_carrier_density_unit_factor::T
    charge_carrier_density_top::T
    charge_carrier_density_bot::T
    n_total_contacts::Int
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

    InvertedCoax{T}() where {T<:AbstractFloat} = new{T}()
end



function InvertedCoax{T}(config_file::Dict)::InvertedCoax where T <: AbstractFloat
    config_file["class"] != "InvertedCoax" ? error() : nothing
    ivc::InvertedCoax = InvertedCoax{T}()
    ivc.name = config_file["name"]
    unit_conversion = Dict{String, T}( "nm"=>1e-9,"um"=>1e-6,"mm"=>1e-3,"cm"=>1e-2,"m"=>1.0 )
    ivc.geometry_unit_factor = unit_conversion[config_file["geometry"]["unit"]]
    materials = Dict("HPGe" => :HPGe, "Vacuum" => :Vacuum)
    ivc.material_detector = material_properties[materials[config_file["materials"]["detector"]]]
    ivc.material_environment = material_properties[materials[config_file["materials"]["environment"]]]
    bulk_types = Dict(  "n" => :ntype,
        "n-type" => :ntype,
        "ntype" => :ntype,
        "p-type" => :ptype,
        "ptype" => :ptype,
        "p" => :ptype  )
    ivc.bulk_type = bulk_types[ config_file["type"] ]
    ivc.cyclic = deg2rad(config_file["cyclic"])
    ivc.mirror_symmetry_θ = config_file["mirror_symmetry_θ"] == "true"
    ivc.borehole_modulation=false
    ivc.crystal_length = ivc.geometry_unit_factor * config_file["geometry"]["crystal"]["length"]
    ivc.crystal_radius = ivc.geometry_unit_factor * config_file["geometry"]["crystal"]["radius"]
    ivc.pc_depth = ivc.geometry_unit_factor * config_file["geometry"]["point_contact"]["depth"]
    ivc.pc_radius = ivc.geometry_unit_factor * config_file["geometry"]["point_contact"]["radius"]
    ivc.groove_depth = ivc.geometry_unit_factor * config_file["geometry"]["groove"]["depth"]
    ivc.groove_rInner = ivc.geometry_unit_factor * config_file["geometry"]["groove"]["rInner"]
    ivc.groove_width = ivc.geometry_unit_factor * config_file["geometry"]["groove"]["width"]
    ivc.groove_endplate = config_file["geometry"]["groove"]["endplate"]
    ivc.borehole_length = ivc.geometry_unit_factor * config_file["geometry"]["borehole"]["length"]
    ivc.borehole_radius = ivc.geometry_unit_factor * config_file["geometry"]["borehole"]["radius"]
    ivc.taper_outer_length = ivc.geometry_unit_factor * config_file["geometry"]["taper"]["outer"]["length"]
    ivc.taper_outer_rInner = ivc.geometry_unit_factor * config_file["geometry"]["taper"]["outer"]["rInner"]
    ivc.taper_outer_angle = atan((ivc.crystal_radius-ivc.taper_outer_rInner)/ivc.taper_outer_length)
    ivc.taper_inner_length = ivc.geometry_unit_factor * config_file["geometry"]["taper"]["inner"]["length"]
    ivc.taper_inner_rOuter = ivc.geometry_unit_factor * config_file["geometry"]["taper"]["inner"]["rOuter"]
    ivc.taper_inner_angle = atan((ivc.taper_inner_rOuter-ivc.borehole_radius)/ivc.taper_inner_length)
    ivc.charge_carrier_density_top = config_file["charge_carrier_density"]["top"]
    ivc.charge_carrier_density_bot = config_file["charge_carrier_density"]["bot"]
    ivc.n_total_contacts = config_file["segmentation"]["n_contacts_total"]
    ivc.n_individual_segments = config_file["segmentation"]["n_individual_segments"]
    ivc.custom_grouping = config_file["segmentation"]["custom_grouping"]

    ivc.volumes=[]
    push!(ivc.volumes,Tubs("Outer_Shape", 1, 16.0, "floating", 0.0,
    0.0,
    ivc.geometry_unit_factor * config_file["geometry"]["crystal"]["radius"],
    0.0,
    2π,
    0.0,
    ivc.geometry_unit_factor * config_file["geometry"]["crystal"]["length"]
    ))

    push!(ivc.volumes,Tubs("Borehole", 2, 1.0, "floating", 0.0,
    0.0,
    ivc.geometry_unit_factor * config_file["geometry"]["borehole"]["radius"],
    0.0,
    2π,
    ivc.geometry_unit_factor * config_file["geometry"]["crystal"]["length"]-ivc.geometry_unit_factor * config_file["geometry"]["borehole"]["length"],
    ivc.geometry_unit_factor * config_file["geometry"]["crystal"]["length"]
    ))

    push!(ivc.volumes,Tubs("Groove", 2, 1.0, "floating", 0.0,
    ivc.geometry_unit_factor * config_file["geometry"]["groove"]["rInner"],
    ivc.geometry_unit_factor * config_file["geometry"]["groove"]["rInner"]+ivc.geometry_unit_factor * config_file["geometry"]["groove"]["width"],
    0.0,
    2π,
    0.0,
    ivc.geometry_unit_factor * config_file["geometry"]["groove"]["depth"]
    ))

    ivc.n_individual_segments >0 ? construct_segmentation_arrays_for_individual_segments(ivc, config_file) : nothing
    construct_floating_boundary_arrays(ivc)
    construct_grouped_channels_array(ivc, config_file)
    return ivc
end

function InvertedCoax(mytype::Type{<:AbstractFloat},config_file::Dict)
    return InvertedCoax{mytype}(config_file)
end

function InvertedCoax{T}(inputfilename::String) where T <: AbstractFloat
    dicttext = read(inputfilename, String)
    parsed_json_file::Dict = JSON.parse(dicttext)
    return InvertedCoax{T}(parsed_json_file)
end

function InvertedCoax(inputfilename::String)
    dicttext = read(inputfilename, String)
    parsed_json_file::Dict = JSON.parse(dicttext)
    return InvertedCoax{Float32}(parsed_json_file)
end

function InvertedCoax(mytype::Type{<:AbstractFloat},inputfilename::String)
    dicttext = read(inputfilename, String)
    parsed_json_file::Dict = JSON.parse(dicttext)
    return InvertedCoax{mytype}(parsed_json_file)
end

