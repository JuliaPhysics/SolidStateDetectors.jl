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

    bulk_geometry::AbstractGeometry{T}

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


function CGD{T}(config_file::Dict)::CGD where {T <: AbstractFloat}
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

    det.bulk_geometry = CartesianBox3D{T}( config_file["geometry"]["crystal"], det.geometry_unit )



    det.charge_carrier_density_top = geom_round(T(config_file["charge_carrier_density"]["top"]))
    det.charge_carrier_density_bot = geom_round(T(config_file["charge_carrier_density"]["bot"]))
    
    det.n_contacts = Contact{T, :N}[ Contact{T, :N}( contact_dict, det.geometry_unit ) for contact_dict in config_file["segmentation"]["n_contacts"]["contacts"] ]
    det.p_contacts = Contact{T, :P}[ Contact{T, :P}( contact_dict, det.geometry_unit ) for contact_dict in config_file["segmentation"]["p_contacts"]["contacts"] ]
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
    println("    crystal: x: $( uconvert.(c.geometry_unit, (c.bulk_geometry.x .* u"m"))) ")
    println("    crystal: y: $( uconvert.(c.geometry_unit, (c.bulk_geometry.y .* u"m"))) ")
    println("    crystal: z: $( uconvert.(c.geometry_unit, (c.bulk_geometry.z .* u"m"))) ")
end
function display(c::CGD)  println(c)   end
function print(c::CGD)  println(c)   end
function show(c::CGD)  println(c)   end
function show(io::IO, ::MIME"text/plain",  c::CGD) 
    show(io, c)
end


function Grid(  det::CGD{T}; 
                init_grid_spacing::Vector{<:Real} = [0.001, 0.001, 0.001], 
                for_weighting_potential::Bool = false)::CartesianGrid3D{T} where {T}

    important_x_points::Vector{T} = T[] #uniq(sort(round.(get_important_r_points(detector), sigdigits=6)))
    important_y_points::Vector{T} = T[] #!only_2d ? sort(get_important_θ_points(detector)) : T[]
    important_z_points::Vector{T} = T[] #uniq(sort(round.(get_important_z_points(detector), sigdigits=6))) #T[]

    init_grid_spacing::Vector{T} = T.(init_grid_spacing)
    
    int_x::Interval{:closed, :closed, T} = Interval{:closed, :closed, T}(det.bulk_geometry.x[1] - 0.005, det.bulk_geometry.x[2] + 0.005)
    int_y::Interval{:closed, :closed, T} = Interval{:closed, :closed, T}(det.bulk_geometry.y[1] - 0.005, det.bulk_geometry.y[2] + 0.005)
    int_z::Interval{:closed, :closed, T} = Interval{:closed, :closed, T}(det.bulk_geometry.z[1] - 0.005, det.bulk_geometry.z[2] + 0.005)
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
    return in(pt, det.bulk_geometry)
end

function get_charge_density(detector::SolidStateDetector{T}, pt::StaticVector{3, T})::T where {T}
    top_net_charge_carrier_density::T = detector.charge_carrier_density_top * 1e10 * 1e6  #  1/cm^3 -> 1/m^3
    bot_net_charge_carrier_density::T = detector.charge_carrier_density_bot * 1e10 * 1e6  #  1/cm^3 -> 1/m^3
    crystal_length_x::T = detector.bulk_geometry.x[2] - detector.bulk_geometry.x[1]
    slope::T = (top_net_charge_carrier_density - bot_net_charge_carrier_density) / crystal_length_x
    ρ::T = bot_net_charge_carrier_density + pt[1] * slope
    return ρ 
end

function get_boundary_value(det::CGD{T}, pt::StaticVector{3, T}, weighting_potential_contact_id::Union{Missing, Int} = missing)::Tuple{Bool, T} where {T}
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
    if is_boundary_point && !ismissing(weighting_potential_contact_id) 
        if contact_id == weighting_potential_contact_id
            boundary_potential = 1
        else
            boundary_potential = 0
        end
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
        grid::Grid{T, 3, :Cartesian}, ssd::CGD{T}; weighting_potential_contact_id::Union{Missing, Int} = missing)::Nothing where {T <: AbstractFloat, N}
    
    is_weighting_potential::Bool = !ismissing(weighting_potential_contact_id)

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
                is_boundary_point, boundary_potential = get_boundary_value(ssd, pt, weighting_potential_contact_id)
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
