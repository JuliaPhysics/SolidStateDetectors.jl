mutable struct InvertedCoax{T<:AbstractFloat} <: SolidStateDetector{T}
    name::String
    material_detector::NamedTuple
    material_environment::NamedTuple
    cyclic::T
    mirror_symmetry_φ::Bool
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
    ivc.mirror_symmetry_φ = config_file["mirror_symmetry_φ"] == "true"
    ivc.borehole_modulation=false
    ivc.crystal_length = geom_round(T(ivc.geometry_unit_factor * config_file["geometry"]["crystal"]["length"]))
    ivc.crystal_radius = geom_round(T(ivc.geometry_unit_factor * config_file["geometry"]["crystal"]["radius"]))
    ivc.pc_depth = geom_round(T(ivc.geometry_unit_factor * config_file["geometry"]["point_contact"]["depth"]))
    ivc.pc_radius = geom_round(T(ivc.geometry_unit_factor * config_file["geometry"]["point_contact"]["radius"]))
    ivc.groove_depth = geom_round(T(ivc.geometry_unit_factor * config_file["geometry"]["groove"]["depth"]))
    ivc.groove_rInner = geom_round(T(ivc.geometry_unit_factor * config_file["geometry"]["groove"]["rInner"]))
    ivc.groove_width = geom_round(T(ivc.geometry_unit_factor * config_file["geometry"]["groove"]["width"]))
    ivc.groove_endplate = config_file["geometry"]["groove"]["endplate"]
    ivc.borehole_length = geom_round(T(ivc.geometry_unit_factor * config_file["geometry"]["borehole"]["length"]))
    ivc.borehole_radius = geom_round(T(ivc.geometry_unit_factor * config_file["geometry"]["borehole"]["radius"]))
    ivc.taper_outer_length = geom_round(T(ivc.geometry_unit_factor * config_file["geometry"]["taper"]["outer"]["length"]))
    ivc.taper_outer_rInner = geom_round(T(ivc.geometry_unit_factor * config_file["geometry"]["taper"]["outer"]["rInner"]))
    ivc.taper_outer_angle = geom_round(T(atan((ivc.crystal_radius-ivc.taper_outer_rInner)/ivc.taper_outer_length)))
    ivc.taper_inner_length = geom_round(T(ivc.geometry_unit_factor * config_file["geometry"]["taper"]["inner"]["length"]))
    ivc.taper_inner_rOuter = geom_round(T(ivc.geometry_unit_factor * config_file["geometry"]["taper"]["inner"]["rOuter"]))
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


function Grid(  detector::Union{Coax{T}, BEGe{T}, InvertedCoax{T}};
                init_grid_spacing::Vector{<:Real} = [0.005, 5.0, 0.005],
                for_weighting_potential::Bool = false)::CylindricalGrid{T} where {T}

    important_r_points::Vector{T} = uniq(sort(round.(get_important_r_points(detector), sigdigits=6)))
    important_φ_points::Vector{T} = T[]#!only_2d ? sort(get_important_φ_points(detector)) : T[]
    important_z_points::Vector{T} = uniq(sort(round.(get_important_z_points(detector), sigdigits=6))) #T[]

    init_grid_spacing::Vector{T} = T.(init_grid_spacing)

    # r
    int_r = Interval{:closed, :closed, T}(0, detector.crystal_radius + 3 * init_grid_spacing[1])
    ax_r::DiscreteAxis{T, :r0, :infinite} = DiscreteAxis{:r0, :infinite}(int_r, step = init_grid_spacing[1])
    rticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_r, important_r_points, atol=init_grid_spacing[1]/4)
    ax_r = DiscreteAxis{T, :r0, :infinite}(int_r, rticks)

    # φ
    if !(0 <= detector.cyclic <= T(2π))
        @warn "'cyclic=$(round(rad2deg(detector.cyclic), digits = 0))°' is ∉ [0, 360] -> 'cyclic' set to 360°."
        detector.cyclic = T(2π)
    end
    int_φ = if !for_weighting_potential
        detector.mirror_symmetry_φ || detector.cyclic == 0 ? Interval{:closed, :closed, T}(0, detector.cyclic / 2) : Interval{:closed, :open, T}(0, detector.cyclic)
    else
        detector.cyclic == 0 ? Interval{:closed, :closed, T}(0, 0) : Interval{:closed, :open, T}(0, 2π)
    end
    nφ::Int = div(int_φ.right, deg2rad(init_grid_spacing[2])) # nφ must be even or equals to 1 ( 1 -> 2D special case)
    if detector.cyclic > 0
        try
            nsym_φ::Int = Int(round(T(2π) / detector.cyclic, digits = 3))
        catch err
            error("360° divided by 'cyclic=$(round(rad2deg(detector.cyclic), digits = 0))°' does not give an integer -> Set 'cyclic' to 360°.")
        end
    end
    ax_φ = if nφ >= 2
        if detector.mirror_symmetry_φ && !for_weighting_potential
            DiscreteAxis{:reflecting, :reflecting}(int_φ, length = iseven(nφ) ? nφ : nφ + 1)
        else
            DiscreteAxis{:periodic, :periodic}(int_φ, length = iseven(nφ) ? nφ : nφ + 1)
        end
    else
        DiscreteAxis{T, :reflecting, :reflecting}(int_φ, T[0])
    end
    if length(ax_φ) > 1
        φticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_φ, important_φ_points, atol=deg2rad(init_grid_spacing[2])/4)
        ax_φ = typeof(ax_φ)(int_φ, φticks)
    end

    #z
    int_z = Interval{:closed, :closed, T}( -3 * init_grid_spacing[3], detector.crystal_length + 3 * init_grid_spacing[3])
    ax_z = DiscreteAxis{:infinite, :infinite}(int_z, step = init_grid_spacing[3])
    zticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_z, important_z_points, atol=init_grid_spacing[3]/2)
    ax_z = typeof(ax_z)(int_z, zticks)
    if isodd(length(ax_z)) # must be even
        int_z = ax_z.interval
        zticks = ax_z.ticks
        push!(zticks, geom_round((zticks[end] + zticks[end-1]) * 0.5))
        sort!(zticks)
        ax_z = DiscreteAxis{T, :infinite, :infinite}(int_z, zticks) # must be even
    end
    @assert iseven(length(ax_z)) "CylindricalGrid must have even number of points in z."

    return CylindricalGrid{T}( (ax_r, ax_φ, ax_z) )
end


function println(io::IO, d::Union{Coax{T}, BEGe{T}, InvertedCoax{T}}) where {T <: AbstractFloat}
    println("________"*d.name*"________\n")
    println("---General Properties---")
    println("Detector Material: \t $(d.material_detector.name)")
    println("Environment Material: \t $(d.material_environment.name)")
    println("Bulk type: \t\t $(d.bulk_type)")
    println("Core Bias Voltage: \t $(d.segment_bias_voltages[1]) V")
    println("Mantle Bias Voltage: \t $(d.segment_bias_voltages[2]) V\n")
    println("---Geometry---")
    println("Outer Crystal Dimensions: ")
    println("Crystal length: \t $(round(d.crystal_length * 1000, sigdigits=6)) mm")
    println("Crystal diameter: \t $(round(2 * d.crystal_radius * 1000, sigdigits=6)) mm")
end

function show(io::IO, d::Union{Coax{T}, BEGe{T}, InvertedCoax{T}}) where {T <: AbstractFloat} println(d) end
function print(io::IO, d::Union{Coax{T}, BEGe{T}, InvertedCoax{T}}) where {T <: AbstractFloat} println(d) end
function display(io::IO, d::Union{Coax{T}, BEGe{T}, InvertedCoax{T}} ) where {T <: AbstractFloat} println(d) end
function show(io::IO,::MIME"text/plain", d::Union{Coax{T}, BEGe{T}, InvertedCoax{T}}) where {T <: AbstractFloat}
    show(io, d)
end
