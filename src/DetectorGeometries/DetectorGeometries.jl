abstract type SolidStateDetector{T} end

abstract type Volume end

struct Tubs{T<:AbstractFloat}<:Volume
    name::String
    hierarchy::Int
    ϵ
    behaviour::String
    potential
    rStart::T
    rStop::T
    θStart::T
    θStop::T
    zStart::T
    zStop::T
    function Tubs(name::String, hierarchy::Int, ϵ, behaviour::String, potential, rStart::T, rStop::T, θStart::T, θStop::T, zStart::T, zStop::T) where T<:AbstractFloat
    return new{T}(name, hierarchy, ϵ, behaviour, potential, rStart, rStop, θStart, θStop, zStart, zStop)
    end
end

function Tubs(name, hierarchy, ϵ, behaviour, potential, rStart::T, rStop::T, θStart::T, θStop::T, zStart::T, zStop::T) where T<:AbstractFloat
    return Tubs{typeof(zStop)}(name, hierarchy, ϵ, behaviour, potential, rStart, rStop, θStart, θStop, zStart, zStop)
end


function uniq(v::Vector{T})::Vector{T} where {T <: Real}
    v1::Vector{T} = Vector{T}()
    if length(v) > 0
        laste::T = v[1]
        push!(v1, laste)
        for e in v
            if e != laste
                laste = e
                push!(v1, laste)
            end
        end
    end
    return v1
end

function merge_axis_ticks_with_important_ticks(ax::DiscreteAxis{T}, impticks::Vector{T}; atol::Real = 0.0001 )::Vector{T} where {T}
    v::Vector{T} = T[]
    for r in impticks push!(v, r) end
    for r in ax push!(v, r) end
    sort!(v)
    v = uniq(v)
    delete_idcs::Vector{Int} = Int[]
    for i in 1:(length(v) - 1)
        if (v[i + 1] - v[i]) < atol
            if !in(v[i], impticks) push!(delete_idcs, i) end
            if !in(v[i + 1], impticks) push!(delete_idcs, i + 1) end
        end
    end
    delete_idcs = sort(uniq(delete_idcs))
    deleteat!(v, delete_idcs) 
    for impv in impticks
        if !in(impv, v)
            error("Important ticks were removed.")
        end
    end
    return v
end

function Grid(  detector::SolidStateDetector{T}; 
                init_grid_spacing::Vector{<:Real} = [0.005, 5.0, 0.005], 
                for_weighting_potential::Bool = false)::CylindricalGrid{T} where {T}

    important_r_points::Vector{T} = uniq(sort(round.(get_important_r_points(detector), sigdigits=6)))
    important_θ_points::Vector{T} = T[]#!only_2d ? sort(get_important_θ_points(detector)) : T[]
    important_z_points::Vector{T} = uniq(sort(round.(get_important_z_points(detector), sigdigits=6))) #T[]

    init_grid_spacing::Vector{T} = T.(init_grid_spacing)
    
    # r
    int_r = Interval{:closed, :closed, T}(0, detector.crystal_radius + 3 * init_grid_spacing[1])
    ax_r::DiscreteAxis{T, :r0, :infinite} = DiscreteAxis{:r0, :infinite}(int_r, step = init_grid_spacing[1])
    rticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_r, important_r_points, atol=init_grid_spacing[1]/4)
    ax_r = DiscreteAxis{T, :r0, :infinite}(int_r, rticks)
    
    # θ
    if !(0 <= detector.cyclic <= T(2π)) 
        @warn "'cyclic=$(round(rad2deg(detector.cyclic), digits = 0))°' is ∉ [0, 360] -> Set 'cyclic' to 360°."
        detector.cyclic = T(2π)
    end
    int_θ = if !for_weighting_potential
        detector.mirror_symmetry_θ > 0 ? Interval{:closed, :closed, T}(0, detector.cyclic / 2) : Interval{:closed, :open, T}(0, detector.cyclic)
    else
        detector.cyclic == 0 ? Interval{:closed, :closed, T}(0, 0) : Interval{:closed, :open, T}(0, 2π)
    end
    nθ::Int = div(int_θ.right, deg2rad(init_grid_spacing[2])) # nθ must be even or equals to 1 ( 1 -> 2D special case)
    if detector.cyclic > 0
        try
            nsym_θ::Int = Int(round(T(2π) / detector.cyclic, digits = 3))
        catch err
            error("360° divided by 'cyclic=$(round(rad2deg(detector.cyclic), digits = 0))°' does not give an integer -> Set 'cyclic' to 360°.")
        end
    end
    ax_θ = nθ >= 2 ? DiscreteAxis{:periodic, :periodic}(int_θ, length = iseven(nθ) ? nθ : nθ + 1) : DiscreteAxis{T, :periodic, :periodic}(int_θ, T[0])
    if length(ax_θ) > 1 
        θticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_θ, important_θ_points, atol=deg2rad(init_grid_spacing[2])/4)
        ax_θ = typeof(ax_θ)(int_θ, θticks)
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

    return CylindricalGrid{T}( (ax_r, ax_θ, ax_z) )
end

include("BEGe.jl")
include("Coax.jl")
include("Inverted_Coax.jl")
# include("DetectorGeometries_V2.jl")

function println(io::IO, d::SolidStateDetector)
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

function show(io::IO, d::SolidStateDetector) println(d) end
function print(io::IO, d::SolidStateDetector) println(d) end
function display(io::IO, d::SolidStateDetector ) println(d) end
function show(io::IO,::MIME"text/plain", d::SolidStateDetector)
    show(io, d)
end

function SolidStateDetector(mytype::Type{<:AbstractFloat},filename::AbstractString)::SolidStateDetector
    SolidStateDetector{mytype}(filename)
end
function SolidStateDetector(filename::AbstractString)::SolidStateDetector
    SolidStateDetector{Float32}(filename)
end

"""
    SolidStateDetector{T}(filename::AbstractString)::SolidStateDetector where T <: AbstractFloat

Reads in a config-JSON file and returns an Detector struct which holds all information specified in the config file.
"""
function SolidStateDetector{T}(filename::AbstractString)::SolidStateDetector where T <: AbstractFloat
    dicttext = read(filename, String)
    parsed_json_file = JSON.parse(dicttext)
    detector_class = parsed_json_file["class"]
    if detector_class == "Coax"
        return Coax{T}(filename)
    elseif detector_class == "BEGe"
        return BEGe{T}(filename)
    elseif detector_class == "InvertedCoax"
        return InvertedCoax{T}(filename)
    else
        error("Config File does not suit any of the predefined detector geometries. You may want to implement your own 'class'")
    end
end

get_precision_type(d::SolidStateDetector{T}) where {T} = T 

function get_important_r_points_from_geometry(c::Coax)
    important_r_points_from_geometry::Vector = []
    push!(important_r_points_from_geometry,c.crystal_radius)
    push!(important_r_points_from_geometry,c.borehole_radius)
    push!(important_r_points_from_geometry,c.taper_inner_bot_rOuter)
    push!(important_r_points_from_geometry,c.taper_inner_top_rOuter)
    important_r_points_from_geometry
end

function get_important_r_points_from_geometry(b::BEGe)
    important_r_points_from_geometry::Vector = []
    push!(important_r_points_from_geometry,b.crystal_radius)
    push!(important_r_points_from_geometry,b.taper_bot_rInner)
    push!(important_r_points_from_geometry,b.taper_top_rInner)
    push!(important_r_points_from_geometry,b.groove_rInner)
    push!(important_r_points_from_geometry,b.groove_rInner+b.groove_width)
    important_r_points_from_geometry
end
function get_important_r_points_from_geometry(ivc::InvertedCoax)
    important_r_points_from_geometry::Vector = []
    for v in ivc.volumes
            push!(important_r_points_from_geometry,v.rStart)
            push!(important_r_points_from_geometry,v.rStop)
    end
    important_r_points_from_geometry
end

function get_important_r_points(d::SolidStateDetector)
    T=get_precision_type(d)
    important_r_points::Vector{T} = []
    ## Outer Shape
    push!(important_r_points,get_important_r_points_from_geometry(d)...)
    ## From Segmentation
    for tuple in d.segmentation_r_ranges
        !in(tuple[1],important_r_points) ? push!(important_r_points,tuple[1]) : nothing
        !in(tuple[2],important_r_points) ? push!(important_r_points,tuple[2]) : nothing
    end
    important_r_points
end

function get_important_θ_points(d::SolidStateDetector)
    T = get_precision_type(d)
    important_θ_points::Vector{T} = []
    for tuple in d.segmentation_phi_ranges
        !in(tuple[1],important_θ_points) ? push!(important_θ_points,tuple[1]) : nothing
        !in(tuple[2],important_θ_points) ? push!(important_θ_points,tuple[2]) : nothing
    end
    for boundary_midpoint in d.segmentation_boundaryMidpoints_horizontal
        !in(boundary_midpoint,important_θ_points) ? push!(important_θ_points,boundary_midpoint) : nothing
    end
    important_θ_points
end

function get_important_z_points_from_geometry(c::Coax)
    important_z_points_from_geometry::Vector = []
    push!(important_z_points_from_geometry,0.0)
    push!(important_z_points_from_geometry,c.crystal_length)
    push!(important_z_points_from_geometry,c.taper_inner_bot_length)
    push!(important_z_points_from_geometry,c.crystal_length-c.taper_inner_top_length)
    important_z_points_from_geometry
end
function get_important_z_points_from_geometry(b::BEGe)
    important_z_points_from_geometry::Vector = []
    push!(important_z_points_from_geometry,0.0)
    push!(important_z_points_from_geometry,b.crystal_length)
    push!(important_z_points_from_geometry,b.taper_bot_length)
    push!(important_z_points_from_geometry,b.crystal_length-b.taper_top_length)
    b.groove_endplate == "top" ? push!(important_z_points_from_geometry,b.crystal_length-b.groove_depth) : push!(important_z_points_from_geometry,b.groove_depth)
    important_z_points_from_geometry
end
function get_important_z_points_from_geometry(ivc::InvertedCoax)
    important_z_points_from_geometry::Vector = []
    for v in ivc.volumes
            push!(important_z_points_from_geometry,v.zStart)
            push!(important_z_points_from_geometry,v.zStop)
    end
    important_z_points_from_geometry
end
function get_important_z_points(d::SolidStateDetector)
    T=get_precision_type(d)
    important_z_points::Vector{T} = []
    ## Outer Shape
    push!(important_z_points,get_important_z_points_from_geometry(d)...)
    ## From Segmentation
    for tuple in d.segmentation_z_ranges
        !in(tuple[1],important_z_points) ? push!(important_z_points,tuple[1]) : nothing
        !in(tuple[2],important_z_points) ? push!(important_z_points,tuple[2]) : nothing
    end
    for boundary_midpoint in d.segmentation_boundaryMidpoints_vertical
        !in(boundary_midpoint,important_z_points) ? push!(important_z_points,boundary_midpoint) : nothing
    end
    important_z_points
end

function construct_segmentation_arrays_from_repetitive_segment(d::SolidStateDetector,config_file::Dict)::Nothing
    n_total_segments::Int = d.n_total_contacts
    T = get_precision_type(d)
    f = d.geometry_unit_factor
    segmentation_r_ranges::Array{Tuple{T,T},1}= []
    segmentation_phi_ranges::Array{Tuple{T,T},1} = []
    segmentation_z_ranges::Array{Tuple{T,T},1} = []
    segment_bias_voltages::Array{T,1} = []
    segmentation_boundaryWidths_horizontal::Array{Tuple{T,T},1} = []
    segmentation_boundaryWidths_vertical::Array{Tuple{T,T},1} = []
    segmentation_boundaryMidpoints_vertical::Array{T,1} = []
    segmentation_boundaryMidpoints_horizontal::Array{T,1} = []

    #Core
    push!(segmentation_r_ranges,(round(config_file["segmentation"]["core"]["rStart"]*f, digits=5),round(config_file["segmentation"]["core"]["rStop"]*f, digits=5)))
    push!(segmentation_phi_ranges,(round(deg2rad(config_file["segmentation"]["core"]["phiStart"]), digits=5),round(deg2rad(config_file["segmentation"]["core"]["phiStop"]), digits=5)))
    push!(segmentation_z_ranges,(round(config_file["segmentation"]["core"]["zStart"]*f, digits=5),round(config_file["segmentation"]["core"]["zStop"]*f, digits=5)))
    push!(segment_bias_voltages,round(config_file["segmentation"]["core"]["potential"], digits=5))

    #repetitive Segments
    # if d.n_repetitive_segments > 0
    rStart=round(config_file["segmentation"]["repetitive_segment"]["rStart"]*f, digits=5)
    rStop=round(config_file["segmentation"]["repetitive_segment"]["rStop"]*f, digits=5)
    boundaryWidth_horizontal = round(deg2rad(config_file["segmentation"]["repetitive_segment"]["boundaryWidth"]["horizontal"]*f*180/(π*d.crystal_radius)), digits=5)
    phiStart=round(deg2rad(config_file["segmentation"]["repetitive_segment"]["phiStart"])+boundaryWidth_horizontal/2, digits=5)
    phiStop=round(deg2rad(config_file["segmentation"]["repetitive_segment"]["phiStop"]) - boundaryWidth_horizontal/2, digits=5)
    zStart=round(config_file["segmentation"]["repetitive_segment"]["zStart"]*f, digits=5)
    boundaryWidth_vertical=round(config_file["segmentation"]["repetitive_segment"]["boundaryWidth"]["vertical"]*f, digits=5)
    zStop=round(config_file["segmentation"]["repetitive_segment"]["zStop"]*f - boundaryWidth_vertical, digits=5)
    potential=round(config_file["segmentation"]["repetitive_segment"]["potential"], digits=5)
    n_vertical_repetitions = config_file["segmentation"]["repetitive_segment"]["repetitions"]["vertical"]
    n_horizontal_repetitions = config_file["segmentation"]["repetitive_segment"]["repetitions"]["horizontal"]
    for v_rseg in 0 : n_vertical_repetitions-1
            for h_rseg in 0 : n_horizontal_repetitions-1
                    push!(segmentation_r_ranges, (rStart,rStop))
                    h_offset::T = h_rseg*(phiStop-phiStart + boundaryWidth_horizontal)
                    push!(segmentation_phi_ranges,(phiStart+h_offset,phiStop+h_offset))
                    push!(segmentation_boundaryMidpoints_horizontal,phiStart+h_offset-boundaryWidth_horizontal/2)
                    v_offset::T=v_rseg*(zStop-zStart + boundaryWidth_vertical)
                    # v_rseg == n_vertical_repetitions -1 ? push!(segmentation_z_ranges,((zStart+v_offset),(zStop+v_offset))) : push!(segmentation_z_ranges,((zStart+v_offset),(zStop+v_offset)))
                    v_rseg == n_vertical_repetitions -1 ? push!(segmentation_z_ranges,((zStart+v_offset),(zStop+v_offset +boundaryWidth_vertical))) : push!(segmentation_z_ranges,((zStart+v_offset),(zStop+v_offset)))
                    push!(segmentation_boundaryMidpoints_vertical,zStart+v_offset + boundaryWidth_vertical/2)
                    push!(segment_bias_voltages,potential)
            end
    end
    d.segmentation_r_ranges   = segmentation_r_ranges
    d.segmentation_phi_ranges = segmentation_phi_ranges
    d.segmentation_z_ranges   = segmentation_z_ranges
    d.segment_bias_voltages = segment_bias_voltages
    d.segmentation_boundaryMidpoints_horizontal = segmentation_boundaryMidpoints_horizontal
    d.segmentation_boundaryMidpoints_vertical = segmentation_boundaryMidpoints_vertical
    nothing
end

function construct_segmentation_arrays_for_individual_segments(d::SolidStateDetector,config_file::Dict)::Nothing
    n_individual_segments::Int = d.n_individual_segments
    T = get_precision_type(d)
    f = d.geometry_unit_factor

    segmentation_r_ranges::Array{Tuple{T,T},1}= []
    segmentation_phi_ranges::Array{Tuple{T,T},1} = []
    segmentation_z_ranges::Array{Tuple{T,T},1} = []
    segment_bias_voltages::Array{T,1} = []
    segmentation_types::Array{String,1} = []
    segmentation_boundaryWidths_horizontal::Array{Tuple{T,T},1} = []
    segmentation_boundaryWidths_vertical::Array{Tuple{T,T},1} = []
    segmentation_boundaryMidpoints_radial::Array{T,1} = []
    segmentation_boundaryMidpoints_vertical::Array{T,1} = []
    segmentation_boundaryMidpoints_horizontal::Array{T,1} = []

    push!(segmentation_r_ranges,(round(config_file["segmentation"]["core"]["rStart"]*f, digits = 5),round(config_file["segmentation"]["core"]["rStop"]*f, digits = 5)))
    push!(segmentation_phi_ranges,(round(deg2rad(config_file["segmentation"]["core"]["phiStart"]), digits = 5),round(deg2rad(config_file["segmentation"]["core"]["phiStop"]), digits = 5)))
    push!(segmentation_z_ranges,(round(config_file["segmentation"]["core"]["zStart"]*f, digits = 5),round(config_file["segmentation"]["core"]["zStop"]*f, digits = 5)))
    push!(segmentation_types,config_file["segmentation"]["core"]["type"])
    push!(segment_bias_voltages,round(config_file["segmentation"]["core"]["potential"], digits = 5))
    for i_idv_seg in 1:n_individual_segments
            ID = "S$i_idv_seg"
            if config_file["segmentation"][ID]["type"] == "Tubs"
            seg_type = "Tubs"
            elseif config_file["segmentation"][ID]["type"] == "Taper"
            seg_type = config_file["segmentation"][ID]["orientation"]
            end
            boundaryWidth_radial = round(config_file["segmentation"][ID]["boundaryWidth"]["radial"]*f, digits = 5)
            rStart=round(config_file["segmentation"][ID]["rStart"]*f, digits = 5)+boundaryWidth_radial
            rStop=round(config_file["segmentation"][ID]["rStop"]*f, digits = 5)
            boundaryWidth_horizontal = round(deg2rad(config_file["segmentation"][ID]["boundaryWidth"]["horizontal"]*f*180/(π*d.crystal_radius)), digits = 5)
            phiStart=round(deg2rad(config_file["segmentation"][ID]["phiStart"]) + boundaryWidth_horizontal/2, digits = 5)
            phiStop=round(deg2rad(config_file["segmentation"][ID]["phiStop"]) - boundaryWidth_horizontal/2, digits = 5)
            boundaryWidth_vertical=round(config_file["segmentation"][ID]["boundaryWidth"]["vertical"]*f, digits = 5)
            zStart=round(config_file["segmentation"][ID]["zStart"]*f + boundaryWidth_vertical/2, digits = 5)
            zStop=round(config_file["segmentation"][ID]["zStop"]*f - boundaryWidth_vertical/2, digits = 5)
            potential=round(config_file["segmentation"][ID]["potential"], digits = 5)
            if config_file["segmentation"][ID]["repetitive"]==true
                n_radial_repetitions = config_file["segmentation"][ID]["repetitions"]["radial"]
                n_vertical_repetitions = config_file["segmentation"][ID]["repetitions"]["vertical"]
                n_horizontal_repetitions = config_file["segmentation"][ID]["repetitions"]["horizontal"]
                for ir in 0:n_radial_repetitions
                    for iv in 0:n_vertical_repetitions
                        for ih in 0:n_horizontal_repetitions
                            r_offset::T = ir*(rStop-rStart + boundaryWidth_radial)
                            h_offset::T = ih*(phiStop-phiStart + boundaryWidth_horizontal)
                            v_offset::T = iv*(zStop-zStart + boundaryWidth_vertical)
                            if iv == 0 && n_vertical_repetitions > 0
                                push!(segmentation_z_ranges,(T(0.0) , zStop + v_offset))

                            elseif iv == n_vertical_repetitions && n_vertical_repetitions > 0
                                push!(segmentation_z_ranges,(zStart + v_offset, d.crystal_length))

                            else
                                push!(segmentation_z_ranges,(zStart + v_offset, zStop + v_offset))
                            end
                            push!(segmentation_r_ranges,(rStart - r_offset, rStop - r_offset))
                            push!(segmentation_phi_ranges,(phiStart + h_offset, phiStop + h_offset))

                            push!(segmentation_types,seg_type)
                            push!(segment_bias_voltages,potential)

                            push!(segmentation_boundaryMidpoints_radial,rStart-r_offset-boundaryWidth_radial/2)
                            push!(segmentation_boundaryMidpoints_horizontal,phiStop+h_offset+boundaryWidth_horizontal/2)

                            iv < n_vertical_repetitions ? push!(segmentation_boundaryMidpoints_vertical,zStop+v_offset+boundaryWidth_vertical/2) : nothing

                        end
                    end
                end
            else
                push!(segmentation_r_ranges, (rStart,rStop))
                push!(segmentation_phi_ranges,(phiStart,phiStop))
                push!(segmentation_z_ranges,(zStart,zStop))
                push!(segmentation_types,seg_type)
                push!(segment_bias_voltages,potential)

                push!(segmentation_boundaryMidpoints_radial,rStart-boundaryWidth_radial/2)
                push!(segmentation_boundaryMidpoints_horizontal,phiStop+boundaryWidth_horizontal/2)
                push!(segmentation_boundaryMidpoints_vertical,zStop+boundaryWidth_vertical/2)
            end

    end
    d.segmentation_r_ranges   = segmentation_r_ranges
    d.segmentation_phi_ranges = segmentation_phi_ranges
    d.segmentation_z_ranges   = segmentation_z_ranges
    d.segmentation_types      = segmentation_types
    d.segment_bias_voltages = segment_bias_voltages
    d.segmentation_boundaryMidpoints_radial = segmentation_boundaryMidpoints_radial
    d.segmentation_boundaryMidpoints_horizontal = segmentation_boundaryMidpoints_horizontal
    d.segmentation_boundaryMidpoints_vertical = segmentation_boundaryMidpoints_vertical
    nothing
end

function construct_floating_boundary_arrays(d::SolidStateDetector)
    ## Groove
    floating_boundary_r_ranges=[]
    floating_boundary_phi_ranges=[]
    floating_boundary_z_ranges=[]
    floating_boundary_types=[]
    push!(floating_boundary_r_ranges,(d.groove_rInner,d.groove_rInner))
    push!(floating_boundary_phi_ranges,(0.0,2π))
    if d.groove_endplate=="bot"
            push!(floating_boundary_z_ranges,(0.0,d.groove_depth))
    else
            push!(floating_boundary_z_ranges,(d.crystal_length-d.groove_depth,d.crystal_length))
    end
    push!(floating_boundary_types,"Tubs")

    push!(floating_boundary_r_ranges,(d.groove_rInner,d.groove_rInner+d.groove_width))
    push!(floating_boundary_phi_ranges,(0.0,2π))
    push!(floating_boundary_z_ranges,(d.groove_depth,d.groove_depth))
    push!(floating_boundary_types,"Tubs")

    push!(floating_boundary_r_ranges,(d.groove_rInner+d.groove_width,d.groove_rInner+d.groove_width))
    push!(floating_boundary_phi_ranges,(0.0,2π))
    if d.groove_endplate=="bot"
            push!(floating_boundary_z_ranges,(0.0,d.groove_depth))
    else
            push!(floating_boundary_z_ranges,(d.crystal_length-d.groove_depth,d.crystal_length))
    end
    push!(floating_boundary_types,"Tubs")


    ## outer_crystal
    push!(floating_boundary_r_ranges,(0.0,d.crystal_radius))
    push!(floating_boundary_phi_ranges,(0.0,2π))
    push!(floating_boundary_z_ranges,(0.0,0.0))
    push!(floating_boundary_types,"Tubs")

    push!(floating_boundary_r_ranges,(d.crystal_radius,d.crystal_radius))
    push!(floating_boundary_phi_ranges,(0.0,2π))
    push!(floating_boundary_z_ranges,(0.0,d.crystal_length))
    push!(floating_boundary_types,"Tubs")

    push!(floating_boundary_r_ranges,(0.0,d.crystal_radius))
    push!(floating_boundary_phi_ranges,(0.0,2π))
    push!(floating_boundary_z_ranges,(d.crystal_length,d.crystal_length))
    push!(floating_boundary_types,"Tubs")

    ## detector specific outer boundaries
    add_detector_specific_outer_boundaries(d,floating_boundary_r_ranges,floating_boundary_phi_ranges,floating_boundary_z_ranges,floating_boundary_types)

    d.floating_boundary_r_ranges=floating_boundary_r_ranges
    d.floating_boundary_phi_ranges=floating_boundary_phi_ranges
    d.floating_boundary_z_ranges=floating_boundary_z_ranges
    d.floating_boundary_types=floating_boundary_types
end

function add_detector_specific_outer_boundaries(c::Coax,rs,phis,zs,mytypes)
    ## taper_top
    if !iszero(c.taper_outer_top_length)
            push!(rs,(c.taper_outer_top_rInner,c.crystal_radius))
            push!(phis,(0.,2π))
            push!(zs,(c.crystal_length-c.taper_outer_top_length,c.crystal_length))
            push!(mytypes,"c//")
    end
    if !iszero(c.taper_inner_top_length)
            push!(rs,(c.taper_inner_top_rOuter,c.crystal_radius))
            push!(phis,(0.,2π))
            push!(zs,(c.crystal_length-c.taper_inner_top_length,c.crystal_length))
            push!(mytypes,"/c")
    end
    ## taper_bot
    if !iszero(c.taper_outer_bot_length)
            push!(rs,(c.taper_outer_bot_rInner,c.crystal_radius))
            push!(phis,(0.,2π))
            push!(zs,(c.taper_outer_bot_length,c.crystal_length))
            push!(mytypes,"c/")
    end
    if !iszero(c.taper_inner_bot_length)
            push!(rs,(c.taper_inner_bot_rOuter,c.crystal_radius))
            push!(phis,(0.,2π))
            push!(zs,(c.taper_inner_bot_length,c.crystal_length))
            push!(mytypes,"//c")
    end
end

function add_detector_specific_outer_boundaries(b::BEGe,rs,phis,zs,mytypes)
    ## taper_top
    if !iszero(b.taper_top_length)
            push!(rs,(b.taper_top_rInner,b.crystal_radius))
            push!(phis,(0.,2π))
            push!(zs,(b.crystal_length-b.taper_top_length,b.crystal_length))
            push!(mytypes,"c//")
    end
    ## taper_bot
    if !iszero(b.taper_bot_length)
            push!(rs,(b.taper_bot_rInner,b.crystal_radius))
            push!(phis,(0.,2π))
            push!(zs,(b.taper_bot_length,b.crystal_length))
            push!(mytypes,"c/")
    end
end
function add_detector_specific_outer_boundaries(ivc::InvertedCoax,rs,phis,zs,mytypes)
        nothing
end


@inline Base.in(point::CylindricalPoint, detector::SolidStateDetector) =
    contains(detector, point)

@inline Base.in(point::StaticArray{Tuple{3},<:Real}, detector::SolidStateDetector) =
    convert(CylindricalPoint, point) in detector

@inline Base.in(point::StaticArray{Tuple{3},<:Quantity}, detector::SolidStateDetector) =
    to_internal_units(u"m", point) in detector

@inline Base.in(point::CoordinateTransformations.Cylindrical, detector::SolidStateDetector) =
    convert(CylindricalPoint, point) in detector


# TODO: Deprecate function contains in favour of Base.in (see above):

# false -> outside
function contains(c::Coax, r::Real, θ::Real, z::Real)::Bool
    rv::Bool = true
    if !check_outer_limits(c, r, θ, z) rv = false end
    if !check_borehole(c, r, θ, z) rv = false end
    if !check_tapers(c, r, θ, z) rv = false end
    return rv
end
function contains(c::Coax, p::CylindricalPoint)::Bool
    rv::Bool = true
    if !check_outer_limits(c, p.r,p.θ,p.z) rv = false end
    if !check_borehole(c, p.r,p.θ,p.z) rv = false end
    if !check_tapers(c, p.r,p.θ,p.z) rv = false end
    return rv
end
function contains(b::BEGe, r::Real, θ::Real, z::Real)::Bool
    rv::Bool = true
    check_outer_limits(b,r,θ,z) ? nothing : rv = false
    check_tapers(b,r,θ,z) ? nothing : rv = false
    check_grooves(b,r,θ,z) ? nothing : rv = false
    rv
end
function contains(b::BEGe, p::CylindricalPoint)::Bool
    rv::Bool = true
    check_outer_limits(b,p.r,p.θ,p.z) ? nothing : rv = false
    check_tapers(b,p.r,p.θ,p.z) ? nothing : rv = false
    check_grooves(b,p.r,p.θ,p.z) ? nothing : rv = false
    rv
end
function contains(b::BEGe, p::Cylindrical)::Bool
    rv::Bool = true
    check_outer_limits(b,p.r,p.θ,p.z) ? nothing : rv = false
    check_tapers(b,p.r,p.θ,p.z) ? nothing : rv = false
    check_grooves(b,p.r,p.θ,p.z) ? nothing : rv = false
    rv
end
function contains(ivc::InvertedCoax,r::Real, θ::Real, z::Real)::Bool
    rv::Bool = true
    check_outer_limits(ivc,r,θ,z) ? nothing : rv = false
    check_tapers(ivc,r,θ,z) ? nothing : rv = false
    if ivc.borehole_modulation == true
        check_borehole(ivc,r,θ,z,ivc.borehole_ModulationFunction) ? nothing : rv = false
    else
        check_borehole(ivc,r,θ,z) ? nothing : rv = false
    end
    check_grooves(ivc,r,θ,z) ? nothing : rv = false
    rv
end
function contains(ivc::InvertedCoax,p::Cylindrical)::Bool
    rv::Bool = true
    check_outer_limits(ivc,p.r,p.θ,p.z) ? nothing : rv = false
    check_tapers(ivc,p.r,p.θ,p.z) ? nothing : rv = false
    if ivc.borehole_modulation == true
        check_borehole(ivc,p.r,p.θ,p.z,ivc.borehole_ModulationFunction) ? nothing : rv = false
    else
        check_borehole(ivc,p.r,p.θ,p.z) ? nothing : rv = false
    end
    check_grooves(ivc,p.r,p.θ,p.z) ? nothing : rv = false
    rv
end
function contains(ivc::InvertedCoax,p::CylindricalPoint)::Bool
    rv::Bool = true
    check_outer_limits(ivc,p.r,p.θ,p.z) ? nothing : rv = false
    check_tapers(ivc,p.r,p.θ,p.z) ? nothing : rv = false
    if ivc.borehole_modulation == true
        check_borehole(ivc,p.r,p.θ,p.z,ivc.borehole_ModulationFunction) ? nothing : rv = false
    else
        check_borehole(ivc,p.r,p.θ,p.z) ? nothing : rv = false
    end
    check_grooves(ivc,p.r,p.θ,p.z) ? nothing : rv = false
    rv
end
# function contains(ivc::InvertedCoax,r::Real, θ::Real, z::Real, accepted_ϵ=[16.0])::Bool
#     rv::Bool = true
#     crystal_basic_shape = ivc.volumes[1]
#     if check_volume(crystal_basic_shape,r,θ,z)
#         for v in ivc.volumes
#             if v.ϵ in accepted_ϵ
#                 nothing
#             else
#                 if check_volume(v,r,θ,z) rv=false end
#             end
#         end
#         check_tapers(ivc,r,θ,z) ? nothing : rv = false
#     else
#         rv = false
#     end
#     rv
# end

function contains(b::BEGe, p::Tuple)::UInt8
    check_outer_limits(b,p) ? nothing : return 0
    check_tapers(b,p) ? nothing : return 0
    check_grooves(b,p) ? nothing : return 0
    return 1
end

function contains(c::Coax, p::Tuple)::UInt8
    check_outer_limits(c,p) ? nothing : return 0
    check_borehole(c,p) ? nothing : return 0
    check_tapers(c,p) ? nothing : return 0
    return 1
end

function check_outer_limits(d::SolidStateDetector, r::Real, θ::Real, z::Real)::Bool
    rv::Bool = true
    if (r > d.crystal_radius) || (z < 0) || (z > d.crystal_length)
        return false
    end
    return rv
end

function check_outer_limits(b::SolidStateDetector, p::Tuple)::Bool
    rv::Bool = true
    T::Type = get_precision_type(b)
    @fastmath begin
        r::T = round(sqrt(p[1]^2 + p[2]^2), digits=5)
        if r > b.crystal_radius || p[3] < 0 || p[3] > b.crystal_length
                rv = false
        end
    end
    return rv
end

function check_borehole(c::Coax, r::Real, θ::Real, z::Real)::Bool
    rv = true
    if r < c.borehole_radius
        rv = false
    end
    return rv
end
function check_borehole(c::Coax,p::Tuple)::Bool
    T::Type = get_precision_type(c)
    rv::Bool = true
    @fastmath begin
        r::T = round(sqrt(p[1]^2 + p[2]^2), digits=5)
        if r < c.borehole_radius
                rv = false
        end
    end
    return rv
end

function check_borehole(ivc::InvertedCoax, r::Real, θ::Real, z::Real)::Bool#returns true if point is not inside borehole
    rv = true
    if r < round(ivc.borehole_radius,digits=6) && z >round(ivc.crystal_length-ivc.borehole_length,digits=6)
        rv = false
    end
    rv
end

function check_borehole(ivc::InvertedCoax, r::T, θ::T, z::T,modulation_function::Function)::Bool where T<:Real#returns true if point is not inside borehole
    rv = true
    epsilon::T=0.000
    # if r < round(ivc.borehole_radius+modulation_function(θ)-epsilon,digits=6) && z >round(ivc.crystal_length-ivc.borehole_length,digits=6)
    if r < round(ivc.borehole_radius+modulation_function(θ),digits=6) && z >round(ivc.crystal_length-ivc.borehole_length,digits=6)
        rv = false
    end
    rv
end

function check_tapers(b::BEGe, p::Tuple)
    T = get_precision_type(b)
    # p[1]=T(p[1])
    # p[2]=T(p[2])
    # p[3]=T(p[3])
    z::T = p[3]
    if z > (b.crystal_length-b.taper_top_length) && z <= b.crystal_length ##Check top taper
        angle_taper_top::T = atan((b.crystal_radius-b.taper_top_rInner)/b.taper_top_length)
        r_taper_top::T = tan(angle_taper_top)*(b.crystal_length-p[3])+b.taper_top_rInner
        r_point::T = sqrt(p[1]^2 + p[2]^2)
        if r_point > r_taper_top
            return false
        else
            nothing
        end
    elseif p[3] < (b.taper_bot_length) && p[3] >= 0.0 ## Check bot taper
        angle_taper_bot = atan((b.crystal_radius-b.taper_bot_rInner)/b.taper_bot_length)
        r_taper_bot = tan(angle_taper_bot) * p[3] + b.taper_bot_rInner
        r_point = rfromxy(p[1],p[2])
        if r_point > r_taper_bot
            # println("bot taper")
            return false
        else
            nothing
        end
    end
    return true
end

function check_tapers(b::BEGe,r::Real,θ::Real,z::Real)::Bool
    rv::Bool = true
    if z > (b.crystal_length-b.taper_top_length) && z <= b.crystal_length ##Check top taper
        angle_taper_top = atan((b.crystal_radius-b.taper_top_rInner)/b.taper_top_length)
        r_taper_top = tan(angle_taper_top)*(b.crystal_length-z)+b.taper_top_rInner
        if r > r_taper_top
            rv = false
        else
            nothing
        end
    elseif z < (b.taper_bot_length) && z >= 0.0 ## Check bot taper
        angle_taper_bot = atan((b.crystal_radius-b.taper_bot_rInner)/b.taper_bot_length)
        r_taper_bot = tan(angle_taper_bot) * z + b.taper_bot_rInner
        if r > r_taper_bot
            rv = false
        else
            nothing
        end
    end
    rv
end

function check_tapers(ivc::InvertedCoax, r::Real, θ::Real, z::Real)::Bool
    rv::Bool = true
    if !iszero(ivc.taper_outer_length) && z > (ivc.crystal_length-ivc.taper_outer_length)
        if r > ivc.crystal_radius-get_r_from_z_for_taper(ivc.taper_outer_angle, z-(ivc.crystal_length-ivc.taper_outer_length))
                rv=false
        end
    end
    if !iszero(ivc.taper_inner_length) && z > (ivc.crystal_length-ivc.taper_inner_length)
        if r < ivc.borehole_radius + get_r_from_z_for_taper(ivc.taper_inner_angle, z-(ivc.crystal_length-ivc.taper_inner_length))
                rv=false
        end
    end
    rv
end

function get_r_from_z_for_taper(angle,z)
    return z*tan(angle)
end

function check_tapers(c::Coax, r::Real, θ::Real, z::Real)::Bool
    if z > (c.crystal_length-c.taper_outer_top_length) && z <= c.crystal_length ##Check top taper
        angle_taper_outer_top = atan((c.crystal_radius-c.taper_outer_top_rInner)/c.taper_outer_top_length)
        r_taper_outer_top = tan(angle_taper_outer_top)*(c.crystal_length-z)+c.taper_outer_top_rInner
        # if r>c.type_precision(r_taper_outer_top)
        if r > r_taper_outer_top
            # println("top outer taper")
            return false
        end
        # elseif z < (c.taper_outer_bot_length) && z >= c.type_precision(0.0) ## Check bot taper
    elseif z < (c.taper_outer_bot_length) && z >= 0.0 ## Check bot taper
        angle_taper_outer_bot = atan((c.crystal_radius-c.taper_outer_bot_rInner)/c.taper_outer_bot_length)
        r_taper_outer_bot = tan(angle_taper_outer_bot) * z + c.taper_outer_bot_rInner
        # if r > c.type_precision(r_taper_outer_bot)
        if r > r_taper_outer_bot
            # println("bot outer taper")
            return false
        end
    end

    if z > (c.crystal_length-c.taper_inner_top_length) && z <= c.crystal_length ##Check top taper
        angle_taper_inner_top = atan((c.taper_inner_top_rOuter-c.borehole_radius)/c.taper_inner_top_length)
        r_taper_inner_top = c.taper_inner_top_rOuter - tan(angle_taper_inner_top) * (c.crystal_length - z)
        # if r < signif(c.type_precision(r_taper_inner_top), 5)
        if r < round(r_taper_inner_top, sigdigits=5)
            # println("top inner taper")
            return false
        end
        # elseif z < (c.taper_inner_bot_length) && z >= c.type_precision(0.0) ## Check bot taper
    elseif z < (c.taper_inner_bot_length) && z >= 0 ## Check bot taper
        angle_taper_inner_bot = atan((c.taper_inner_bot_rOuter-c.borehole_radius)/c.taper_inner_bot_length)
        r_taper_inner_bot = c.taper_inner_bot_rOuter - tan(angle_taper_inner_bot) * z
        # if r < signif(c.type_precision(r_taper_inner_bot), 5)
        if r < round(r_taper_inner_bot, sigdigits=5)
            # println("bot inner taper")
            return false
        end
    end
    return true
end

function check_tapers(c::Coax, p::Tuple)::Bool
        T::Type = get_precision_type(c)
        rv::Bool = true
        @fastmath begin
        if p[3] > (c.crystal_length - c.taper_outer_top_length) && p[3] <= c.crystal_length ##Check top taper
            angle_taper_outer_top::T = atan((c.crystal_radius - c.taper_outer_top_rInner) / c.taper_outer_top_length)
            r_taper_outer_top::T = tan(angle_taper_outer_top) * (c.crystal_length - p[3]) + c.taper_outer_top_rInner
            r_point::T = round(sqrt(p[1]^2 + p[2]^2), digits=5)
            if r_point > r_taper_outer_top
                # println("top outer taper")
                rv = false
            end
        elseif p[3] < c.taper_outer_bot_length && p[3] >= 0 ## Check bot taper
            angle_taper_outer_bot = atan((c.crystal_radius - c.taper_outer_bot_rInner) / c.taper_outer_bot_length)
            r_taper_outer_bot = tan(angle_taper_outer_bot) * p[3] + c.taper_outer_bot_rInner
            r_point = round(sqrt(p[1]^2 + p[2]^2), digits=5)
            if r_point > r_taper_outer_bot
                # println("bot outer taper")
                rv = false
            end
        end

        if p[3] > (c.crystal_length - c.taper_inner_top_length) && p[3] <= c.crystal_length ##Check top taper
            angle_taper_inner_top = atan((c.taper_inner_top_rOuter - c.borehole_radius) / c.taper_inner_top_length)
            r_taper_inner_top = c.taper_inner_top_rOuter - tan(angle_taper_inner_top) * (c.crystal_length - p[3])
            r_point = round(sqrt(p[1]^2 + p[2]^2), digits=5)
            if r_point < r_taper_inner_top
                # println("top inner taper")
                rv = false
            end
        elseif p[3] < (c.taper_inner_bot_length) && p[3] >= T(0) ## Check bot taper
            angle_taper_inner_bot = atan((c.taper_inner_bot_rOuter - c.borehole_radius) / c.taper_inner_bot_length)
            r_taper_inner_bot = c.taper_inner_bot_rOuter - tan(angle_taper_inner_bot) * p[3]
            r_point = round(sqrt(p[1]^2 + p[2]^2), digits=5)
            if r_point < r_taper_inner_bot
                # println("bot inner taper")
                rv = false
            end
        end
end
return rv
end

function check_grooves(b::BEGe,r::Real,θ::Real,z::Real)::Bool
    # T = get_precision_type(b)
    rv::Bool=true
    if b.groove_endplate == "bot"
        if z >= 0 && z < b.groove_depth
            if r > b.groove_rInner && r < (b.groove_rInner+b.groove_width)
                rv = false
            end
        end
    elseif b.groove_endplate == "top"
        if z < b.crystal_length && z > (b.crystal_length-b.groove_depth)
            if r > b.groove_rInner && r < (b.groove_rInner+b.groove_width)
                rv = false
            end
        end
    end
    rv
end

function check_grooves(ivc::InvertedCoax,r::Real,θ::Real,z::Real)::Bool
    T=typeof(z)
    rv::Bool=true
    if z >= T(0.0) && z < ivc.groove_depth
        if r>ivc.groove_rInner && r < (ivc.groove_rInner+ivc.groove_width)
            rv = false
        else
            nothing
        end
    end
    rv
end

function check_volume(v::Volume,r::Real,θ::Real,z::Real)::Bool #return true if in volume
    rv=true
    T=typeof(v.rStart)
    if typeof(v) == Tubs{T}
        if r< v.rStart || r > v.rStop
            rv = false
        elseif θ < v.θStart || θ > v.θStop
            rv = false
        elseif z < v.zStart || z > v.zStop
            rv = false
        end
    end
    rv
end

function check_grooves(b::BEGe, p::Tuple)
    T = get_precision_type(b)
    if p[3] > T(0) && p[3]<b.groove_depth
        r_point::T = rfromxy(p[1],p[2])
        if r_point>b.groove_rInner && r_point< (b.groove_rInner+b.groove_width)
            # println("groove")
            return false
        else
            nothing
        end
    end
    return true
end

function rfromxy(x::Real, y::Real)
    return sqrt(x^2+y^2)
end

# is_boundary
# function is_boundary_point(d::SolidStateDetector, r::Real, θ::Real, z::Real)::Bool
#     rv::Bool = false
#     for iseg in 1:size(d.segment_bias_voltages,1)
#         if (d.segmentation_r_ranges[iseg][1] <= r <= d.segmentation_r_ranges[iseg][2])
#             if (d.segmentation_phi_ranges[iseg][1] <= θ <= d.segmentation_phi_ranges[iseg][2])
#                 if (d.segmentation_z_ranges[iseg][1] <= z <= d.segmentation_z_ranges[iseg][2])
#                     rv = true
#                 end
#             end
#         end
#     end
#     rv
# end
# function is_boundary_point(d::Coax, r::Real, θ::Real, z::Real, rs, θs, zs)
#     is_boundary_point(d,r,θ,z)
# end
function is_boundary_point(d::SolidStateDetector, r::T, θ::T, z::T, rs::Vector{T}, θs::Vector{T}, zs::Vector{T}) where T <:AbstractFloat
    rv::Bool = false
    if θ < 0
        θ += d.cyclic
    end
    digits::Int=6
    if  d.borehole_modulation == true
        idx_r_closest_gridpoint_to_borehole = searchsortednearest(rs, d.borehole_radius + d.borehole_ModulationFunction(θ))
        idx_r = findfirst(x->x==r,rs)
        bore_hole_r = rs[idx_r_closest_gridpoint_to_borehole]
    else
        nothing
    end
    
    for iseg in 1:size(d.segment_bias_voltages,1)
        if d.borehole_modulation == true && iseg == d.borehole_segment_idx
             # x = (round(d.segmentation_r_ranges[iseg][1]+d.borehole_ModulationFunction(θ)-ϵ,digits=digits) <= r <= round(d.segmentation_r_ranges[iseg][2]+d.borehole_ModulationFunction(θ)+ϵ,digits=digits))
             x = (idx_r_closest_gridpoint_to_borehole == idx_r)
        elseif d.borehole_modulation == true && iseg == d.borehole_bot_segment_idx
             x = (d.segmentation_r_ranges[iseg][1] <= r <= round(bore_hole_r,digits=digits))
             # x = idx_r_closest_gridpoint_to_borehole == idx_r
        elseif d.borehole_modulation == true && iseg == d.borehole_top_segment_idx
            x = (round(bore_hole_r,digits=digits) <= r <= d.segmentation_r_ranges[iseg][2])
            # x = idx_r_closest_gridpoint_to_borehole == idx_r
        else
            x = (d.segmentation_r_ranges[iseg][1] <= r <= d.segmentation_r_ranges[iseg][2])
        end
        if x
            if (d.segmentation_phi_ranges[iseg][1] <= θ <= d.segmentation_phi_ranges[iseg][2])
                if (d.segmentation_z_ranges[iseg][1] <= z <= d.segmentation_z_ranges[iseg][2])
                    if d.segmentation_types[iseg]=="Tubs"
                        rv = true
                    else
                        if isapprox(rs[searchsortednearest(rs,analytical_taper_r_from_θz(θ,z,
                            d.segmentation_types[iseg],
                            d.segmentation_r_ranges[iseg],
                            d.segmentation_phi_ranges[iseg],
                            d.segmentation_z_ranges[iseg]
                            ))],r)
                            rv = true
                        end
                    end
                end
            end
        end
    end
    rv
end

function point_type(d::SolidStateDetector, p::Cylindrical)
    point_type(d,p.r,p.θ,p.z)
end
function point_type(d::SolidStateDetector, r::T, θ::T, z::T) where {T<:AbstractFloat}
    T==Float32 ? atol = 0.000001 : atol = 0.000000000000001
    rv::Symbol = :bulk
    !contains(d, Cylindrical{T}(r,θ,z)) ? rv = :outside : nothing
    i=0
    digits::Int=6
    ############################# Electrode Definitions
    for iseg in 1:size(d.segment_bias_voltages,1)
        if d.borehole_modulation == true && iseg == d.borehole_segment_idx
             x = (round(d.segmentation_r_ranges[iseg][1]+d.borehole_ModulationFunction(θ),digits=digits) <= r <= round(d.segmentation_r_ranges[iseg][2]+d.borehole_ModulationFunction(θ),digits=digits))
        elseif d.borehole_modulation == true && iseg == d.borehole_bot_segment_idx
             x = (d.segmentation_r_ranges[iseg][1] <= r <= round(d.segmentation_r_ranges[iseg][2]+d.borehole_ModulationFunction(θ),digits=digits))
        elseif d.borehole_modulation == true && iseg == d.borehole_top_segment_idx
            x = (round(d.segmentation_r_ranges[iseg][1]+d.borehole_ModulationFunction(θ),digits=digits) <= r <= d.segmentation_r_ranges[iseg][2])
        else
            x = (d.segmentation_r_ranges[iseg][1] <= r <= d.segmentation_r_ranges[iseg][2])
        end
        if x
            if (d.segmentation_phi_ranges[iseg][1] <= θ <= d.segmentation_phi_ranges[iseg][2])
                if (d.segmentation_z_ranges[iseg][1] <= z <= d.segmentation_z_ranges[iseg][2])
                    if d.segmentation_types[iseg]=="Tubs"
                        rv = :electrode
                        i=iseg
                    else
                        if isapprox(analytical_taper_r_from_θz(θ,z,
                            d.segmentation_types[iseg],
                            d.segmentation_r_ranges[iseg],
                            d.segmentation_phi_ranges[iseg],
                            d.segmentation_z_ranges[iseg]
                            ),r, atol = atol )
                        rv = :electrode
                        i=iseg
                        end
                    end
                end
            end
        end
    end

    if i==0
    ############################ Floating Boundary Definitions
    for iseg in 1:size(d.floating_boundary_r_ranges,1)
        if (d.floating_boundary_r_ranges[iseg][1] <= r <= d.floating_boundary_r_ranges[iseg][2])
            if (d.floating_boundary_phi_ranges[iseg][1] <= θ <= d.floating_boundary_phi_ranges[iseg][2])
                if (d.floating_boundary_z_ranges[iseg][1] <= z <= d.floating_boundary_z_ranges[iseg][2])
                    if d.floating_boundary_types[iseg]=="Tubs"
                        rv = :floating_boundary
                        i=iseg
                    else
                        if isapprox(analytical_taper_r_from_θz(θ,z,
                            d.floating_boundary_types[iseg],
                            d.floating_boundary_r_ranges[iseg],
                            d.floating_boundary_phi_ranges[iseg],
                            d.floating_boundary_z_ranges[iseg]
                            ),r, atol = atol )
                        rv = :floating_boundary
                        i=iseg
                    end
                end
            end
        end
    end
end
end
rv, i
end

function analytical_taper_r_from_θz(θ,z,orientation,r_bounds,θ_bounds,z_bounds)
    r=0
    angle = atan(abs(r_bounds[2]-r_bounds[1])/abs(z_bounds[2]-z_bounds[1]))
    if orientation == "c//"
        r = r_bounds[2] - tan(angle)*(z-minimum(z_bounds))
    elseif orientation == "/c"
        r = r_bounds[1] + tan(angle)*(z-minimum(z_bounds))
    end
    r
end



function get_segment_idx(d::SolidStateDetector,r::Real,θ::Real,z::Real)
    digits=6
    for iseg in 1:size(d.segment_bias_voltages,1)
        if d.borehole_modulation == true && iseg == d.borehole_segment_idx
             x = (round(d.segmentation_r_ranges[iseg][1]+d.borehole_ModulationFunction(θ),digits=digits) <= r <= round(d.segmentation_r_ranges[iseg][2]+d.borehole_ModulationFunction(θ),digits=digits))
        elseif d.borehole_modulation == true && iseg == d.borehole_bot_segment_idx
             x = (d.segmentation_r_ranges[iseg][1] <= r <= round(d.segmentation_r_ranges[iseg][2]+d.borehole_ModulationFunction(θ),digits=digits))
        elseif d.borehole_modulation == true && iseg == d.borehole_top_segment_idx
            x = (round(d.segmentation_r_ranges[iseg][1]+d.borehole_ModulationFunction(θ),digits=digits) <= r <= d.segmentation_r_ranges[iseg][2])
        else
            x = (d.segmentation_r_ranges[iseg][1] <= r <= d.segmentation_r_ranges[iseg][2])
        end
        if x && (d.segmentation_phi_ranges[iseg][1] <= θ <= d.segmentation_phi_ranges[iseg][2]) && (d.segmentation_z_ranges[iseg][1] <= z <= d.segmentation_z_ranges[iseg][2])
            return iseg
        else
            nothing
        end
    end
    return -1
end
function get_segment_idx(d::SolidStateDetector,r::Real,θ::Real,z::Real,rs::Vector{<:Real})
    digits=6
    idx_r_closest_gridpoint_to_borehole = searchsortednearest(rs,d.borehole_radius+d.borehole_ModulationFunction(θ))
    idx_r = findfirst(x->x==r,rs)
    bore_hole_r = rs[idx_r_closest_gridpoint_to_borehole]
    for iseg in 1:size(d.segment_bias_voltages,1)
        if d.borehole_modulation == true && iseg == d.borehole_segment_idx
             # x = (round(d.segmentation_r_ranges[iseg][1]+d.borehole_ModulationFunction(θ)-ϵ,digits=digits) <= r <= round(d.segmentation_r_ranges[iseg][2]+d.borehole_ModulationFunction(θ)+ϵ,digits=digits))
             x = (idx_r_closest_gridpoint_to_borehole == idx_r)
        elseif d.borehole_modulation == true && iseg == d.borehole_bot_segment_idx
             x = (d.segmentation_r_ranges[iseg][1] <= r <= round(bore_hole_r,digits=digits))
             # x = idx_r_closest_gridpoint_to_borehole == idx_r
        elseif d.borehole_modulation == true && iseg == d.borehole_top_segment_idx
            x = (round(bore_hole_r,digits=digits) <= r <= d.segmentation_r_ranges[iseg][2])
            # x = idx_r_closest_gridpoint_to_borehole == idx_r
        else
            x = (d.segmentation_r_ranges[iseg][1] <= r <= d.segmentation_r_ranges[iseg][2])
        end
        if x && (d.segmentation_phi_ranges[iseg][1] <= θ <= d.segmentation_phi_ranges[iseg][2]) && (d.segmentation_z_ranges[iseg][1] <= z <= d.segmentation_z_ranges[iseg][2])
            return iseg
        else
            nothing
        end
    end
    return -1
end

# get_potential
function get_boundary_value(d::SolidStateDetector, r::Real, θ::Real, z::Real, rs::Vector{<:Real})::get_precision_type(d)
    if d.borehole_modulation==true
        res::get_precision_type(d)=0.0
        try
            res = d.segment_bias_voltages[ get_segment_idx(d, r, θ, z, rs) ]
        catch
            println("kek")
            @show r, θ, z
            res = d.segment_bias_voltages[end]
        end
        return res
    else
        return d.segment_bias_voltages[ get_segment_idx(d, r, θ, z) ]
    end
end



function get_charge_density(detector::SolidStateDetector, r::Real, θ::Real, z::Real)::get_precision_type(detector)
    T = get_precision_type(detector)
    top_net_charge_carrier_density::T = detector.charge_carrier_density_top * 1e10 * 1e6  #  1/cm^3 -> 1/m^3
    bot_net_charge_carrier_density::T = detector.charge_carrier_density_bot * 1e10 * 1e6  #  1/cm^3 -> 1/m^3
    slope::T = (top_net_charge_carrier_density - bot_net_charge_carrier_density) / detector.crystal_length
    ρ::T = bot_net_charge_carrier_density + z * slope
    # ρ = ρ - (r - detector.borehole_radius) * 0.2 * ρ / (detector.crystal_radius - detector.borehole_radius) 
    return ρ 
end

function json_to_dict(inputfile::String)::Dict
    parsed_json_file = Dict()
    open(inputfile,"r") do f
        global parsed_json_file
        dicttext = readstring(f)
        parsed_json_file = JSON.parse(dicttext)
    end
    return parsed_json_file
end
include("plot_recipes.jl")
