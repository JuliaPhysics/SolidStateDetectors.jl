mutable struct BEGe{T<:AbstractFloat} <: SolidStateDetector{T}
#### Importet Values from JSON file
    name::String
    material_detector::NamedTuple
    material_environment::NamedTuple
    cyclic::T
    mirror_symmetry_θ::Bool
    bulk_type::Symbol
    geometry_unit_factor::T
    crystal_length::T
    crystal_radius::T
    pc_depth::T
    pc_radius::T
    groove_endplate::String
    groove_depth::T
    groove_rInner::T
    groove_width::T
    taper_bot_length::T
    taper_bot_rInner::T
    taper_top_length::T
    taper_top_rInner::T
    sfc_depth::T
    borehole_modulation::Bool

    charge_carrier_density_unit_factor::T
    charge_carrier_density_top::T
    charge_carrier_density_bot::T
    #Segmentation
    n_total_contacts::Int
    # n_repetitive_segments::Int
    n_individual_segments::Int
    n_total_segment_pieces::Int
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
    BEGe{T}() where {T<:AbstractFloat} = new{T}()
end



# function add_point_to_shape(shape::Array{Array{T,1},1},tuple::Tuple{T,T})::Nothing  where T<:AbstractFloat
#     push!(shape[1],tuple[1])
#     push!(shape[2],tuple[2])
#     nothing
# end

function BEGe(mytype::Type,config_file::Dict)
    return BEGe{mytype}(config_file)
end

function BEGe{T}(config_file::Dict)::BEGe where T<:Real
    println(config_file["class"])
    if config_file["class"] != "BEGe"
    error()
    end
    b::BEGe = BEGe{T}()
    b.name = config_file["name"]
    unit_conversion = Dict{String, T}( "nm"=>1e-9,"um"=>1e-6,"mm"=>1e-3,"cm"=>1e-2,"m"=>1.0f0 )
    b.geometry_unit_factor = unit_conversion[config_file["geometry"]["unit"]]
    materials = Dict("HPGe" => :HPGe, "Vacuum" => :Vacuum)
    b.material_detector = material_properties[materials[config_file["materials"]["detector"]]]
    b.material_environment = material_properties[materials[config_file["materials"]["environment"]]]
    bulk_types = Dict(  "n" => :ntype,
        "n-type" => :ntype,
        "ntype" => :ntype,
        "p-type" => :ptype,
        "ptype" => :ptype,
        "p" => :ptype  )
    b.cyclic = deg2rad(config_file["cyclic"])
    b.mirror_symmetry_θ = config_file["mirror_symmetry_θ"] == "true"
    b.bulk_type = bulk_types[ config_file["type"]  ]
    b.borehole_modulation=false
    b.crystal_length = round(b.geometry_unit_factor * config_file["geometry"]["crystal"]["length"], sigdigits=5)
    b.crystal_radius = round(b.geometry_unit_factor * config_file["geometry"]["crystal"]["radius"], sigdigits=5)
    b.pc_depth = round(b.geometry_unit_factor * config_file["geometry"]["point_contact"]["depth"], sigdigits=5)
    b.pc_radius = round(b.geometry_unit_factor * config_file["geometry"]["point_contact"]["radius"], sigdigits=5)
    b.groove_endplate = config_file["geometry"]["groove"]["endplate"]
    b.groove_depth = b.geometry_unit_factor * config_file["geometry"]["groove"]["depth"]
    b.groove_rInner = b.geometry_unit_factor * config_file["geometry"]["groove"]["rInner"]
    b.groove_width = b.geometry_unit_factor * config_file["geometry"]["groove"]["width"]
    b.taper_bot_length = b.geometry_unit_factor * config_file["geometry"]["taper"]["bot"]["length"]
    b.taper_bot_rInner = b.geometry_unit_factor * config_file["geometry"]["taper"]["bot"]["rInner"]
    b.taper_top_length = b.geometry_unit_factor * config_file["geometry"]["taper"]["top"]["length"]
    b.taper_top_rInner = b.geometry_unit_factor * config_file["geometry"]["taper"]["top"]["rInner"]
    b.charge_carrier_density_top = round(config_file["charge_carrier_density"]["top"], digits = 5)
    b.charge_carrier_density_bot = round(config_file["charge_carrier_density"]["bot"], digits = 5)
    b.n_total_contacts = config_file["segmentation"]["n_contacts_total"]

    b.n_individual_segments = config_file["segmentation"]["n_individual_segments"]
    b.custom_grouping = config_file["segmentation"]["custom_grouping"]

    b.n_individual_segments >0 ? construct_segmentation_arrays_for_individual_segments(b,config_file) : nothing
    construct_floating_boundary_arrays(b)
    construct_grouped_channels_array(b,config_file)
    return b
end

function BEGe(inputfilename::String)
    dicttext = read(inputfilename, String)
    parsed_json_file = JSON.parse(dicttext)
    return BEGe{Float32}(parsed_json_file)
end
function BEGe(mytype::Type{<:AbstractFloat},inputfilename::String)
    dicttext = read(inputfilename, String)
    parsed_json_file = JSON.parse(dicttext)
    return BEGe{mytype}(parsed_json_file)
end
function BEGe{T}(inputfilename::String) where T <: AbstractFloat
    dicttext = read(inputfilename, String)
    parsed_json_file = JSON.parse(dicttext)
    return BEGe{T}(parsed_json_file)
end
function construct_grouped_channels_array(d,config_file)
    channels::Array{Array{Int,1},1}=[]
    if d.custom_grouping==true
        for ichn in 1:d.n_total_contacts
            parts_belonging_to_channel::String = config_file["segmentation"]["Chn$ichn"]
            push!(channels,map(x->parse(Int,x),split(parts_belonging_to_channel,"-")))
        end
    else
        push!(channels,[1])
        push!(channels,[i for i in 2:d.n_total_contacts])
    end
    d.grouped_channels=channels
end

