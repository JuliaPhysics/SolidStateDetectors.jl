

abstract type AbstractGeometry{T} end
mutable struct CartesianBox3D{T} <: AbstractGeometry{T}
    x::Tuple{T, T}
    y::Tuple{T, T}
    z::Tuple{T, T}
end

function in(pt::StaticVector{3, T}, g::CartesianBox3D{T})::Bool where {T}
    return (g.x[1] <= pt[1] <= g.x[2]) && (g.y[1] <= pt[2] <= g.y[2]) && (g.z[1] <= pt[3] <= g.z[2])
end

abstract type AbstractContact{T} end
mutable struct Contact{T, D} <: AbstractContact{T}
    potential::T
    id::Int
    geometry::AbstractGeometry{T}
end

function in(pt::StaticVector{3, T}, c::Contact{T})::Bool where {T}
    return in(pt, c.geometry)
end

# CGD: Cuboid Germanium Detector
mutable struct CGD{T<:AbstractFloat} <: SolidStateDetector{T}
    name::String
    material_detector::NamedTuple
    material_environment::NamedTuple

    mirror_axes_pos::NTuple{3, T}

    path_to_configfile::String
    #### Importet Values from JSON file
    ### Geometry
    bulk_type::Symbol
    
    geometry_unit::Unitful.Units

    crystal_length_x::T
    crystal_length_y::T
    crystal_length_z::T

    charge_carrier_density_unit_factor::T
    charge_carrier_density_top::T
    charge_carrier_density_bot::T

    ### Segmentaion
    n_contacts::Vector{Contact{T, :N}}
    p_contacts::Vector{Contact{T, :P}}

    CGD{T}() where {T<:AbstractFloat} = new{T}()
end

function CGD(mytype::Type{<:AbstractFloat},config_file::Dict)::CGD
    return  CGD{mytype}(config_file::Dict)
end


function CGD{T}(config_file::Dict)::CGD where T<:AbstractFloat
    config_file["class"] != "CGD" ? error() : nothing
    det::CGD = CGD{T}()
    det.name = config_file["name"]
    unit_conversion = Dict{String, Unitful.Units}( "nm" => u"nm", "um" => u"μm", "mm" => u"mm", "cm" => u"cm", "m" => u"m")
    det.geometry_unit = unit_conversion[config_file["geometry"]["unit"]]
    materials = Dict("HPGe" => :HPGe, "Vacuum" => :Vacuum)
    det.material_detector = material_properties[materials[config_file["materials"]["detector"]]]
    det.material_environment = material_properties[materials[config_file["materials"]["environment"]]]
    bulk_types = Dict{Any, Symbol}( "n" => :ntype,
        "n-type" => :ntype,
        "ntype" => :ntype,
        "p-type" => :ptype,
        "ptype" => :ptype,
        "p" => :ptype  )
    det.bulk_type = bulk_types[ config_file["type"] ]

    crystal_dimensions = config_file["geometry"]["crystal"]["dimensions"]

    det.crystal_length_x = geom_round(ustrip(uconvert(u"m", T(crystal_dimensions["x"]) * det.geometry_unit)))
    det.crystal_length_y = geom_round(ustrip(uconvert(u"m", T(crystal_dimensions["x"]) * det.geometry_unit)))
    det.crystal_length_z = geom_round(ustrip(uconvert(u"m", T(crystal_dimensions["y"]) * det.geometry_unit)))

    det.charge_carrier_density_top = geom_round(T(config_file["charge_carrier_density"]["top"]))
    det.charge_carrier_density_bot = geom_round(T(config_file["charge_carrier_density"]["bot"]))

    n_contacts::Vector{Contact{T, :N}} = Contact{T, :N}[]
    p_contacts::Vector{Contact{T, :P}} = Contact{T, :P}[]

    json_n_contact_array = config_file["segmentation"]["n_contacts"]["contacts"]
    json_p_contact_array = config_file["segmentation"]["p_contacts"]["contacts"]

    for c in json_n_contact_array
        new_contact::Contact{T, :N} = Contact{T, :N}(
            c["potential"],
            c["id"],
            CartesianBox3D{T}(
                (geom_round(ustrip(uconvert(u"m", T(c["xStart"]) * det.geometry_unit ))), geom_round(ustrip(uconvert(u"m", T(c["xStop"]) * det.geometry_unit)))),
                (geom_round(ustrip(uconvert(u"m", T(c["yStart"]) * det.geometry_unit ))), geom_round(ustrip(uconvert(u"m", T(c["yStop"]) * det.geometry_unit)))),
                (geom_round(ustrip(uconvert(u"m", T(c["zStart"]) * det.geometry_unit ))), geom_round(ustrip(uconvert(u"m", T(c["zStop"]) * det.geometry_unit))))
            )
        )
        push!(n_contacts, new_contact)
    end
    for c in json_p_contact_array
        new_contact::Contact{T, :P} = Contact{T, :P}(
            c["potential"],
            c["id"],
            CartesianBox3D{T}(
                (geom_round(ustrip(uconvert(u"m", T(c["xStart"]) * det.geometry_unit))), geom_round(ustrip(uconvert(u"m", T(c["xStop"]) * det.geometry_unit)))),
                (geom_round(ustrip(uconvert(u"m", T(c["yStart"]) * det.geometry_unit))), geom_round(ustrip(uconvert(u"m", T(c["yStop"]) * det.geometry_unit)))),
                (geom_round(ustrip(uconvert(u"m", T(c["zStart"]) * det.geometry_unit))), geom_round(ustrip(uconvert(u"m", T(c["zStop"]) * det.geometry_unit))))
            )
        )
        push!(p_contacts, new_contact)
    end
    det.n_contacts = n_contacts
    det.p_contacts = p_contacts
    return det
end

function CGD(mytype::Type{<:AbstractFloat},inputfilename::String)::CGD
    dicttext = read(inputfilename, String)
    parsed_json_file = JSON.parse(dicttext)
    return CGD{mytype}(parsed_json_file)
end
function CGD{T}(inputfilename::String)::CGD where T <: AbstractFloat
    dicttext = read(inputfilename, String)
    parsed_json_file = JSON.parse(dicttext)
    return CGD{T}(parsed_json_file)
end

function CGD(inputfilename::String)::CGD
    dicttext = read(inputfilename, String)
    parsed_json_file = JSON.parse(dicttext)
    return CGD{Float32}(parsed_json_file)
end

function println(c::CGD)
    println("Cuboid Germanium Detector:")
    println("  Geometry:")
    println("    crystal length x = $( uconvert(c.geometry_unit, (c.crystal_length_x * u"m"))) ")
    println("    crystal length y = $( uconvert(c.geometry_unit, (c.crystal_length_y * u"m"))) ")
    println("    crystal length z = $( uconvert(c.geometry_unit, (c.crystal_length_z * u"m"))) ")
end
function show(c::CGD)  println(c)   end
function display(c::CGD)  println(c)   end
function print(c::CGD)  println(c)   end


function Grid(  det::CGD{T}; 
                init_grid_spacing::Vector{<:Real} = [0.001, 0.001, 0.001], 
                for_weighting_potential::Bool = false)::CartesianGrid3D{T} where {T}

    important_x_points::Vector{T} = T[] #uniq(sort(round.(get_important_r_points(detector), sigdigits=6)))
    important_y_points::Vector{T} = T[] #!only_2d ? sort(get_important_θ_points(detector)) : T[]
    important_z_points::Vector{T} = T[] #uniq(sort(round.(get_important_z_points(detector), sigdigits=6))) #T[]

    init_grid_spacing::Vector{T} = T.(init_grid_spacing)
    
    int_x::Interval{:closed, :closed, T} = Interval{:closed, :closed, T}(-det.crystal_length_x, 2 * det.crystal_length_x)
    int_y::Interval{:closed, :closed, T} = Interval{:closed, :closed, T}(-det.crystal_length_y, 2 * det.crystal_length_y)
    int_z::Interval{:closed, :closed, T} = Interval{:closed, :closed, T}(-det.crystal_length_z, 2 * det.crystal_length_z)
    ax_x::DiscreteAxis{T, :infinite, :infinite} = DiscreteAxis{:infinite, :infinite}(int_z, step = init_grid_spacing[1]) 
    ax_y::DiscreteAxis{T, :infinite, :infinite} = DiscreteAxis{:infinite, :infinite}(int_z, step = init_grid_spacing[2]) 
    ax_z::DiscreteAxis{T, :infinite, :infinite} = DiscreteAxis{:infinite, :infinite}(int_z, step = init_grid_spacing[3]) 

    if isodd(length(ax_x)) # RedBlack dimension must be of even length
        xticks = ax_x.ticks
        push!(xticks, geom_round((xticks[end] + xticks[end-1]) * 0.5))
        sort!(xticks)
        ax_x = DiscreteAxis{T, :infinite, :infinite}(int_x, xticks) # must be even
    end
    @assert iseven(length(ax_x)) "CartesianGrid3D must have even number of points in z."
    
    return CartesianGrid3D{T}( (ax_x, ax_y, ax_z) ) 
end


function in(pt::StaticVector{3, T}, det::CGD{T})::Bool where {T}
    return (0 <= pt[1] <= det.crystal_length_x) && (0 <= pt[2] <= det.crystal_length_y) && (0 <= pt[3] <= det.crystal_length_z)
end

function get_charge_density(detector::SolidStateDetector{T}, pt::StaticVector{3, T})::T where {T}
    top_net_charge_carrier_density::T = detector.charge_carrier_density_top * 1e10 * 1e6  #  1/cm^3 -> 1/m^3
    bot_net_charge_carrier_density::T = detector.charge_carrier_density_bot * 1e10 * 1e6  #  1/cm^3 -> 1/m^3
    slope::T = (top_net_charge_carrier_density - bot_net_charge_carrier_density) / detector.crystal_length_z
    ρ::T = bot_net_charge_carrier_density + pt[3] * slope
    return ρ 
end

function get_boundary_value(det::CGD{T}, pt::StaticVector{3, T}, weighting_potential_channel_ids::Union{Missing, Vector{Int}} = missing)::Tuple{Bool, T} where {T}
    is_boundary_point::Bool = false
    boundary_potential::T = 0
    contact_id::Int = -1
    for contact in det.n_contacts
        if in(pt, contact)
            is_boundary_point = true
            contact_id = contact.id
            boundary_potential = contact.potential
            break
        end
    end
    if !is_boundary_point 
        for contact in det.p_contacts
            if in(pt, contact)
                is_boundary_point = true
                contact_id = contact.id
                boundary_potential = contact.potential
                break
            end
        end
    end
    if !ismissing(weighting_potential_channel_ids) 
        error("todo...")
    end
    return is_boundary_point, boundary_potential
end

function get_ρ_and_ϵ(pt::StaticVector{3, T}, ssd::CGD{T})::Tuple{T, T} where {T <: AbstractFloat}
    if in(pt, ssd)
        ρ::T = get_charge_density(ssd, pt) * elementary_charge
        ϵ::T = ssd.material_detector.ϵ_r
        return ρ, ϵ
    else
        ρ = 0
        ϵ = ssd.material_environment.ϵ_r 
        return ρ, ϵ
    end
end

function set_pointtypes_and_fixed_potentials!(pointtypes::Array{PointType, N}, potential::Array{T, N}, 
        grid::Grid{T, 3, :Cartesian}, ssd::CGD{T}; weighting_potential_channel_idx::Union{Missing, Int} = missing)::Nothing where {T <: AbstractFloat, N}
    
    is_weighting_potential::Bool = !ismissing(weighting_potential_channel_idx)

    weighting_potential_channel_ids = if is_weighting_potential
        error("to be implemented")
        # ssd.grouped_channels[weighting_potential_channel_idx]
    else
        missing
    end

    ax_x::Vector{T} = grid[:x].ticks
    ax_y::Vector{T} = grid[:y].ticks
    ax_z::Vector{T} = grid[:z].ticks
    for iz in axes(potential, 3)
        z::T = ax_z[iz]
        for iy in axes(potential, 2)
            y::T = ax_y[iy]
            for ix in axes(potential, 1)
                x::T = ax_x[ix]
                pt::StaticVector{3, T} = SVector{3, T}( x, y, z )              
                is_boundary_point, boundary_potential = get_boundary_value(ssd, pt, weighting_potential_channel_ids)
                if is_boundary_point                    
                    potential[ ix, iy, iz ]  = boundary_potential
                    pointtypes[ ix, iy, iz ] = zero(PointType)
                elseif in(pt, ssd)
                    pointtypes[ ix, iy, iz ] += pn_junction_bit 
                end

            end
        end
    end
    nothing
end
