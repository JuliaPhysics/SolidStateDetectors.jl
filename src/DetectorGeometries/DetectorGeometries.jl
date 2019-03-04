# """
#     abstract type SolidStateDetector{T} <: AbstractConfig{T}
#
# Supertype of all detector structs.
# """
abstract type CGD{T <: AbstractFloat} <: AbstractConfig{T} end
# abstract type SolidStateDetector{T <: AbstractFloat} <: AbstractConfig{T} end



bulk_types = Dict("n" => :ntype,
    "n-type" => :ntype,
    "ntype" => :ntype,
    "p-type" => :ptype,
    "ptype" => :ptype,
    "p" => :ptype  )
unit_conversion = Dict{String, Unitful.Units}( "nm" => u"nm", "um" => u"μm", "mm" => u"mm", "cm" => u"cm", "m" => u"m")


include("Geometries/Geometries.jl")
include("Contacts.jl")

include("SolidStateDetector.jl")
# Cartesian:
# include("CGD.jl")
get_precision_type(d::SolidStateDetector{T}) where {T} = T

"""
    SolidStateDetector{T}(filename::AbstractString)::SolidStateDetector{T} where {T <: AbstractFloat}

Reads in a config-JSON file and returns an Detector struct which holds all information specified in the config file.
"""
function parse_config_file(filename::AbstractString)::Dict{Any,Any}
    if endswith(filename, ".toml")
        error("currently only .json and .yaml files are supported. We intend to add .toml support in the near future")
    elseif endswith(filename, ".json")
        dicttext = read(filename, String)
        parsed_dict = JSON.parse(dicttext)
    elseif endswith(filename, ".yaml")
        parsed_dict = YAML.load_file(filename)
    else
        error("currently only .json and .yaml files are supported.")
    end
    parsed_dict
end

function yaml2json_convert(filename::String)
    data = YAML.load(open(filename))
    open(replace(filename, ".yaml" => ".json"),"w") do f
            JSON.print(f, data, 4)
    end
end
function yaml2json(directory::String)# or filename
    if endswith(directory, ".yaml")
        yaml2json_convert(directory)
    else
        yamlfiles = filter(x->endswith(x,".yaml"),readdir(directory))
        for filename in yamlfiles
            yaml2json_convert(filename)
        end
    end
end

function SolidStateDetector{T}(filename::AbstractString)::SolidStateDetector{T} where {T <: AbstractFloat}
    parsed_dict = parse_config_file(filename)
    detector_class = parsed_dict["class"]
    global unit_conversions = Dict{String,T}( "nm"=>1e-9,"um"=>1e-6,"mm"=>1e-3,"cm"=>1e-2,"m"=>1.0 )
    if detector_class == "Coax"
        return SolidStateDetector{T}(parsed_dict)
    elseif detector_class == "BEGe"
        return SolidStateDetector{T}(parsed_dict)
    elseif detector_class == "InvertedCoax"
        return SolidStateDetector{T}(parsed_dict)
    elseif detector_class == "CGD"
        return SolidStateDetector{T}(parsed_dict)
    elseif detector_class == "CustomDetector"
        return SolidStateDetector{T}(parsed_dict)
    else
        error("Config File does not suit any of the predefined detector geometries. You may want to implement your own 'class'")
    end
end
function SolidStateDetector(T::Type{<:AbstractFloat} = Float32, filename::AbstractString = SSD_examples[:InvertedCoax])::SolidStateDetector{T}
    SolidStateDetector{T}(filename)
end
function SolidStateDetector(filename::AbstractString)::SolidStateDetector{Float32}
    SolidStateDetector{Float32}(filename)
end

@inline in(point::CylindricalPoint, detector::SolidStateDetector) =
    contains(detector, point)

@inline in(point::StaticArray{Tuple{3}, <:Real}, detector::SolidStateDetector) =
    convert(CylindricalPoint, point) in detector

@inline in(point::StaticArray{Tuple{3}, <:Quantity}, detector::SolidStateDetector) =
    to_internal_units(u"m", point) in detector

@inline in(point::CartesianPoint, detector::SolidStateDetector) =
    convert(CylindricalPoint, point) in detector



function get_important_r_points(c::SolidStateDetector{T}) where T
    r::Vector{T}=[]
    for contact in c.contacts
        for g in contact.geometry
            push!(r, get_r(g)...)
        end
    end
    sort!(unique!(r))
end

function get_important_φ_points(c::SolidStateDetector{T}) where T
    φ::Vector{T}=[]
    for contact in c.contacts
        for g in contact.geometry
            push!(φ, get_φ(g)...)
        end
    end
    sort!(unique!(φ))
end

function get_important_z_points(c::SolidStateDetector{T}) where T
    z::Vector{T}=[]
    for contact in c.contacts
        for g in contact.geometry
            push!(z, get_z(g)...)
        end
    end
    sort!(unique!(z))
end

function is_boundary_point(c::SolidStateDetector, r::T, φ::T, z::T, rs::Vector{T}, φs::Vector{T}, zs::Vector{T})::Tuple{Bool,Real,Int} where T <:AbstractFloat
    p = CylindricalPoint{T}(r, φ, z)
    if false #!(p in c)
        return false, 0.0, 0
    else
        for contact in c.contacts
            if in(p, contact, rs)
                return true, contact.potential, contact.id
            end
        end
        for external_part in c.external_parts
            if in(p, external_part, rs)
                return true, external_part.potential, external_part.id
            end
        end 
        return false, 0.0, 0
    end
end

function point_type(c::SolidStateDetector{T}, p::CylindricalPoint{T})::Tuple{Symbol,Int} where T
    for contact in c.contacts
        if p in contact || geom_round(p) in contact || in(go_to_nearest_gridpoint(c,p), contact, c.rs) || in(go_to_nearest_gridpoint(c,geom_round(p)), contact, c.rs) 
            return :electrode, contact.id
        end
    end
    if is_surface_point(c,p)[1] 
        for contact in c.contacts
            if in(go_to_nearest_gridpoint(c,p), contact, c.rs) 
                return :electrode, contact.id
            else
                return :floating_boundary, 0
            end
        end    
    elseif !(p in c)
        return :outside, -1
    else
        return :bulk, -1
    end
end

function go_to_nearest_gridpoint(d::SolidStateDetector{T}, p::CylindricalPoint{T})::CylindricalPoint{T} where T
    CylindricalPoint{T}(d.rs[searchsortednearest(d.rs,p.r)],d.φs[searchsortednearest(d.φs,p.φ)],d.zs[searchsortednearest(d.zs,p.z)])
end

# function is_surface_point(c::SolidStateDetector{T}, p::CylindricalPoint{T})::Bool where T
#     if !(p in c)
#         return false
#     elseif !(
#         CylindricalPoint{T}(prevfloat(p.r),p.φ,p.z) in c
#         && CylindricalPoint{T}(nextfloat(p.r),p.φ,p.z) in c
#         && CylindricalPoint{T}(p.r,prevfloat(p.φ),p.z) in c
#         && CylindricalPoint{T}(p.r,nextfloat(p.φ),p.z) in c
#         && CylindricalPoint{T}(p.r,p.φ,prevfloat(p.z)) in c
#         && CylindricalPoint{T}(p.r,p.φ,nextfloat(p.z)) in c)
#         return true
#     else
#         return false
#     end
# end

function is_surface_point(c::SolidStateDetector{T}, p::CylindricalPoint{T})::Tuple{Bool,CartesianPoint{T}} where T
    if !(p in c)
        return false, CartesianPoint{T}(0.0,0.0,0.0)
    end
    n::MVector{3,T} = @MVector T[0.0,0.0,0.0]
    look_around::Vector{Bool} = [CylindricalPoint{T}(prevfloat(p.r),p.φ,p.z) in c,
        CylindricalPoint{T}(nextfloat(p.r),p.φ,p.z) in c,
        CylindricalPoint{T}(p.r,prevfloat(p.φ),p.z) in c,
        CylindricalPoint{T}(p.r,nextfloat(p.φ),p.z) in c, 
        CylindricalPoint{T}(p.r,p.φ,prevfloat(p.z)) in c, 
        CylindricalPoint{T}(p.r,p.φ,nextfloat(p.z)) in c]
    if !(false in look_around)
        return false , CartesianPoint{T}(n...)
    else
        look_around[1]==false ? n[1] -= 1 : nothing
        look_around[2]==false ? n[1] += 1 : nothing
        look_around[3]==false ? n[2] -= 1 : nothing
        look_around[4]==false ? n[2] += 1 : nothing
        look_around[5]==false ? n[3] -= 1 : nothing
        look_around[6]==false ? n[3] += 1 : nothing 
        # println(look_around , " " , n)
        Rα::SMatrix{3,3,T} = @SArray([cos(p.φ) -1*sin(p.φ) 0;sin(p.φ) cos(p.φ) 0;0 0 1])
        return true, geom_round(CartesianPoint((Rα * n)...))
    end
end

function get_charge_density(detector::SolidStateDetector{T}, r::Real, φ::Real, z::Real)::T where {T <: AbstractFloat}
    top_net_charge_carrier_density::T = detector.charge_carrier_density_top * 1e10 * 1e6  #  1/cm^3 -> 1/m^3
    bot_net_charge_carrier_density::T = detector.charge_carrier_density_bot * 1e10 * 1e6  #  1/cm^3 -> 1/m^3
    slope::T = (top_net_charge_carrier_density - bot_net_charge_carrier_density) / detector.crystal_length
    ρ::T = bot_net_charge_carrier_density + z * slope
    return ρ
end


function get_ρ_and_ϵ(pt::CylindricalPoint{T}, ssd::SolidStateDetector{T})::Tuple{T, T} where {T <: AbstractFloat}
    if in(pt, ssd)
        ρ::T = get_charge_density(ssd, pt.r, pt.φ, pt.z) * elementary_charge
        ϵ::T = ssd.material_detector.ϵ_r
        return ρ, ϵ
    elseif in(pt,ssd.external_parts)
        for ep in ssd.external_parts
            if pt in ep
                ρ = 0
                ϵ = ep.material.ϵ_r
            end
        end
        return ρ, ϵ
    else
        ρ = 0
        ϵ = ssd.material_environment.ϵ_r
        return ρ, ϵ
    end
end

function write_grid_to_detector(ssd::SolidStateDetector{T}, grid::Grid{T, 3, :Cylindrical}) where {T}
    ssd.rs = grid[:r].ticks 
    ssd.φs = grid[:φ].ticks
    ssd.zs = grid[:z].ticks
end

function set_pointtypes_and_fixed_potentials!(pointtypes::Array{PointType, N}, potential::Array{T, N},
        grid::Grid{T, N, :Cylindrical}, ssd::SolidStateDetector{T}; weighting_potential_contact_id::Union{Missing, Int} = missing)::Nothing where {T <: AbstractFloat, N}

    channels::Array{Int, 1} = if !ismissing(weighting_potential_contact_id)
        [weighting_potential_contact_id]
    else
        Int[]
    end

    axr::Vector{T} = grid[:r].ticks
    axφ::Vector{T} = grid[:φ].ticks
    axz::Vector{T} = grid[:z].ticks

    for iz in axes(potential, 3)
        z::T = axz[iz]
        for iφ in axes(potential, 2)
            φ::T = axφ[iφ]
            for ir in axes(potential, 1)
                r::T = axr[ir]
                pt::CylindricalPoint{T} = CylindricalPoint{T}( r, φ, z )
                b, boundary_potential, contact_id = is_boundary_point(ssd, r, φ, z, axr, axφ, axz)
                if b
                    pot::T = if ismissing(weighting_potential_contact_id)
                        boundary_potential
                    else
                        contact_id == weighting_potential_contact_id ? 1 : 0
                    end
                    potential[ ir, iφ, iz ] = pot
                    pointtypes[ ir, iφ, iz ] = zero(PointType)
                elseif in(pt, ssd)
                    pointtypes[ ir, iφ, iz ] += pn_junction_bit
                end

            end
        end
    end
    nothing
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
function bounding_box(d::SolidStateDetector{T})::NamedTuple where T
    (
    r_range = ClosedInterval{T}(0.0,d.crystal_radius),
    φ_range = ClosedInterval{T}(0.0,2π),
    z_range = ClosedInterval{T}(0.0,d.crystal_length)
    )
end

function Grid(  detector::SolidStateDetector{T};
                init_grid_spacing::Vector{<:Real} = [0.005, 5.0, 0.005],
                for_weighting_potential::Bool = false)::CylindricalGrid{T} where {T}

    important_r_points::Vector{T} = uniq(sort(round.(get_important_r_points(detector), sigdigits=6)))
    important_φ_points::Vector{T} = T[]#!only_2d ? sort(get_important_φ_points(detector)) : T[]
    important_z_points::Vector{T} = uniq(sort(round.(get_important_z_points(detector), sigdigits=6))) #T[]

    init_grid_spacing::Vector{T} = T.(init_grid_spacing)

    # r
    int_r = Interval{:closed, :closed, T}(detector.world.r_interval)
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
    int_z = Interval{:closed, :closed, T}( detector.world.z_interval)
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

include("plot_recipes.jl")
