

abstract type AbstractGeometry{T} end
mutable struct CartesianBox3D{T} <: AbstractGeometry{T}
    x::Tuple{T, T}
    y::Tuple{T, T}
    z::Tuple{T, T}
end

abstract type AbstractContact{T} end
mutable struct Contact{T, D} <: AbstractContact{T}
    potential::T
    id::Int
    geometry::AbstractGeometry{T}
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
    println()
    println("\tCuboid Germanium Detector:")
    println("\t++++++++++++++++++++++++++")
    println("\tGeometry:")
    println("\t  length x = $( uconvert(c.geometry_unit, (c.crystal_length_x * u"m"))) ")
    println("\t  length y = $( uconvert(c.geometry_unit, (c.crystal_length_y * u"m"))) ")
    println("\t  length z = $( uconvert(c.geometry_unit, (c.crystal_length_z * u"m"))) ")
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

    if isodd(length(ax_z)) # RedBlack dimension must be of even length
        zticks = ax_z.ticks
        push!(zticks, geom_round((zticks[end] + zticks[end-1]) * 0.5))
        sort!(zticks)
        ax_z = DiscreteAxis{T, :infinite, :infinite}(int_z, zticks) # must be even
    end
    @assert iseven(length(ax_z)) "CartesianGrid3D must have even number of points in z."
    
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