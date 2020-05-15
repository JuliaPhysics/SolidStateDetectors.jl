"""
    mutable struct SolidStateDetector{T <: SSDFloat, CS} <: AbstractConfig{T}

CS: Coordinate System: -> :cartesian / :cylindrical
"""
mutable struct SolidStateDetector{T <: SSDFloat, CS} <: AbstractConfig{T}
    name::String  # optional
    inputunits::Dict{String, Unitful.Units}
    world::World{T, 3}

    config_dict::Dict

    medium::NamedTuple # this should become a struct at some point

    semiconductors::Vector{Semiconductor{T}}
    contacts::Vector{Contact{T}}
    passives::Vector{Passive{T}}

    virtual_drift_volumes::Vector{AbstractVirtualVolume{T}}
end

get_precision_type(d::SolidStateDetector{T}) where {T} = T
get_coordinate_system(d::SolidStateDetector{T, CS}) where {T, CS} = CS

function SolidStateDetector{T, S}()::SolidStateDetector{T} where {T <: SSDFloat, S}
    semiconductors::Vector{Semiconductor{T}}, contacts::Vector{Contact{T}}, passives::Vector{Passive{T}} = [], [], []
    virtual_drift_volumes::Vector{AbstractVirtualVolume{T}} = []
    world_limits = get_world_limits_from_objects(Val(S), semiconductors, contacts, passives )
    world = World(Val(S), world_limits)

    return SolidStateDetector{T, S}(
        "EmptyDetector",
        default_unit_dict(),
        world,
        Dict(),
        material_properties[materials["vacuum"]],
        semiconductors,
        contacts,
        passives,
        virtual_drift_volumes
    )
end

function SolidStateDetector{T}()::SolidStateDetector{T} where {T <: SSDFloat}
    S::Symbol = :cartesian
    return SolidStateDetector{T, S}()
end
function SolidStateDetector()::SolidStateDetector{Float32, :cartesian}
    return SolidStateDetector{Float32, :cartesian}()
end

function default_unit_dict()::Dict{String, Unitful.Units}
    return Dict{String, Unitful.Units}(
        "length" => u"m", # change this to u"m" ? SI Units
        "potential" => u"V",
        "angle" => u"°",
        "temperature" => u"K"
    )
end


function construct_units(config_file_dict::Dict)
    dunits::Dict{String, Unitful.Units} = default_unit_dict()
    if haskey(config_file_dict, "units")
        d = config_file_dict["units"]
        if haskey(d, "length") dunits["length"] = unit_conversion[d["length"]] end
        if haskey(d, "angle") dunits["angle"] = unit_conversion[d["angle"]] end
        if haskey(d, "potential") dunits["potential"] = unit_conversion[d["potential"]] end
        if haskey(d, "temperature") dunits["temperature"] = unit_conversion[d["temperature"]] end
    end
    dunits
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

function construct_virtual_volume(T, pass::Dict, inputunit_dict::Dict{String, Unitful.Units})
    construct_virtual_volume(T, pass, inputunit_dict, Val{Symbol(pass["model"])} )
end
function construct_virtual_volume(T, pass::Dict, inputunit_dict::Dict{String, Unitful.Units}, ::Type{Val{:dead}})
    DeadVolume{T}(pass, inputunit_dict)
end
function construct_virtual_volume(T, pass::Dict, inputunit_dict::Dict{String, Unitful.Units}, ::Type{Val{:arbitrary}})
    ArbitraryDriftModificationVolume{T}(pass, inputunit_dict)
end

function construct_objects(T, objects::Vector, semiconductors, contacts, passives, virtual_drift_volumes, inputunit_dict)::Nothing
    for obj in objects
        if obj["type"] == "semiconductor"
            push!(semiconductors, construct_semiconductor(T, obj, inputunit_dict))
        elseif obj["type"] == "contact"
            push!(contacts, construct_contact(T, obj, inputunit_dict))
        elseif obj["type"] == "passive"
            push!(passives, construct_passive(T, obj, inputunit_dict))
        elseif obj["type"] == "virtual_drift_volume"
            push!(virtual_drift_volumes, construct_virtual_volume(T, obj, inputunit_dict))
        else
            @warn "please specify the class to be either a \"semiconductor\", a \"contact\", or \"passive\""
        end
    end
    nothing
end

function get_world_limits_from_objects(S::Val{:cylindrical}, s::Vector{Semiconductor{T}}, c::Vector{Contact{T}}, p::Vector{Passive{T}}) where {T <: SSDFloat}
    ax1l::T, ax1r::T, ax2l::T, ax2r::T, ax3l::T, ax3r::T = 0, 1, 0, 1, 0, 1
    imps_1::Vector{T} = []
    imps_3::Vector{T} = []
    for objects in [s, c, p]
        for object in objects
            for posgeo in object.geometry_positive
                append!(imps_1, get_important_points( posgeo, Val{:r}()))
                append!(imps_3, get_important_points( posgeo, Val{:z}()))
            end
        end
    end
    imps_1 = uniq(sort(imps_1))
    imps_3 = uniq(sort(imps_3))
    if length(imps_1) > 1
        ax1l = minimum(imps_1)
        ax1r = maximum(imps_1)
    elseif length(imps_1) == 1
        ax1l = minimum(imps_1)
        ax1r = maximum(imps_1) + 1
    end
    if length(imps_3) > 1
        ax3l = minimum(imps_3)
        ax3r = maximum(imps_3)
    elseif length(imps_3) == 1
        ax3l = minimum(imps_3)
        ax3r = maximum(imps_3) + 1
    end
    return ax1l, ax1r, ax2l, ax2r, ax3l, ax3r
end
function get_world_limits_from_objects(S::Val{:cartesian}, s::Vector{Semiconductor{T}}, c::Vector{Contact{T}}, p::Vector{Passive{T}}) where {T <: SSDFloat}
    ax1l::T, ax1r::T, ax2l::T, ax2r::T, ax3l::T, ax3r::T = 0, 1, 0, 1, 0, 1
    imps_1::Vector{T} = []
    imps_2::Vector{T} = []
    imps_3::Vector{T} = []
    for objects in [s, c, p]
        for object in objects
            for posgeo in object.geometry_positive
                append!(imps_1, get_important_points( posgeo, Val{:x}()))
                append!(imps_2, get_important_points( posgeo, Val{:y}()))
                append!(imps_3, get_important_points( posgeo, Val{:z}()))
            end
        end
    end
    imps_1 = uniq(sort(imps_1))
    imps_2 = uniq(sort(imps_2))
    imps_3 = uniq(sort(imps_3))
    if length(imps_1) > 1
        ax1l = minimum(imps_1)
        ax1r = maximum(imps_1)
    elseif length(imps_1) == 1
        ax1l = minimum(imps_1)
        ax1r = maximum(imps_1) + 1
    end
    if length(imps_2) > 1
        ax2l = minimum(imps_2)
        ax2r = maximum(imps_2)
    elseif length(imps_2) == 1
        ax2l = minimum(imps_2)
        ax2r = maximum(imps_2) + 1
    end
    if length(imps_3) > 1
        ax3l = minimum(imps_3)
        ax3r = maximum(imps_3)
    elseif length(imps_3) == 1
        ax3l = minimum(imps_3)
        ax3r = maximum(imps_3) + 1
    end
    return ax1l, ax1r, ax2l, ax2r, ax3l, ax3r
end

function SolidStateDetector{T}(config_file::Dict)::SolidStateDetector{T} where{T <: SSDFloat}
    grid_type::Symbol = :cartesian
    semiconductors::Vector{Semiconductor{T}}, contacts::Vector{Contact{T}}, passives::Vector{Passive{T}} = [], [], []
    virtual_drift_volumes::Vector{AbstractVirtualVolume{T}} = []
    medium::NamedTuple = material_properties[materials["vacuum"]]
    inputunits = dunits::Dict{String, Unitful.Units} = Dict{String, Unitful.Units}(
        "length" => u"m", # change this to u"m" ? SI Units
        "potential" => u"V",
        "angle" => u"°",
        "temperature" => u"K"
    )
    inputunits = construct_units(config_file)
    if haskey(config_file, "medium")
        medium = material_properties[materials[config_file["medium"]]]
    end

    if haskey(config_file, "objects")
        construct_objects(T, config_file["objects"], semiconductors, contacts, passives, virtual_drift_volumes, inputunits)
    end

    world = if haskey(config_file, "grid")
        if isa(config_file["grid"], Dict)
            grid_type = Symbol(config_file["grid"]["coordinates"])
            World(T, config_file["grid"], inputunits)
        elseif isa(config_file["grid"], String)
            grid_type = Symbol(config_file["grid"])
            world_limits = get_world_limits_from_objects(Val(grid_type), semiconductors, contacts, passives)
            World(Val(grid_type), world_limits)
        end
    else
        world_limits = get_world_limits_from_objects(Val(grid_type), semiconductors, contacts, passives )
        World(Val(grid_type), world_limits)
    end

    c = SolidStateDetector{T, grid_type}()
    c.name = haskey(config_file, "name") ? config_file["name"] : "NoNameDetector"
    c.config_dict = config_file
    c.semiconductors = semiconductors
    c.contacts = contacts
    c.passives = passives
    c.inputunits = inputunits
    c.medium = medium
    c.world = world
    c.virtual_drift_volumes = virtual_drift_volumes
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

function println(io::IO, d::SolidStateDetector{T, CS}) where {T <: SSDFloat, CS}
    println("________"*d.name*"________\n")
    # println("Class: ",d.class)
    println("---General Properties---")
    println("- Precision type: $(T)")
    println("- Environment Material: $(d.medium.name)")
    println("- Grid Type: $(CS)")
    println()
    println("# Semiconductors: $(length(d.semiconductors))")
    for (isc, sc)  in enumerate(d.semiconductors)
        println("\t_____Semiconductor $(isc)_____\n")
        println(sc)
    end
    println()
    println("# Contacts: $(length(d.contacts))")
    if length(d.contacts)<=5
        for c in d.contacts
            println(c)
        end
    end
    println()
    println("# Passives: $(length(d.passives))")
    if length(d.passives)<=5
        for p in d.passives
            # println(c)
        end
    end
end

function show(io::IO, d::SolidStateDetector{T}) where {T <: SSDFloat} println(io, d) end
function print(io::IO, d::SolidStateDetector{T}) where {T <: SSDFloat} println(io, d) end
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



function paint_object(det::SolidStateDetector{T}, object::AbstractObject, grid::Grid{T, 3, S}) where {T <: SSDFloat, S}
    samples = []
    axes_syms = S == :cylindrical ? [:r, :φ, :z] : [:x, :y, :z]
    all_imps = [get_important_points(det, s) for s in axes_syms]
    for g in object.geometry_positive
        stepsizes = T[]
        for (iax, sax) in enumerate(axes_syms)
            imps_ax = all_imps[iax]
            ax::Vector{T} = grid.axes[iax].ticks
            imin::Int = searchsortednearest(ax, minimum(imps_ax))
            imax::Int = searchsortednearest(ax, maximum(imps_ax))
            
            ax = ax[imin:imax]
            delete_inds::Vector{Int} = Int[]
            for imp in imps_ax
                push!(delete_inds, searchsortednearest(ax, imp))
            end
            unique!(sort!(delete_inds))
            deleteat!(ax, delete_inds)
            stepsize::T = length(ax) <= 1 ? T(1) : geom_round((minimum(diff(ax)) / 4))
            imps_g = get_important_points(g, Val(sax))
            unique!(sort!(imps_g))
            if length(imps_g) > 1
                min_imps_g::T = minimum(diff(imps_g)) / 4
                if min_imps_g < stepsize
                    stepsize = geom_round(min_imps_g)
                end
            end
            imps_g_min::T, imps_g_max::T = if length(imps_g) > 0
                minimum(imps_g), maximum(imps_g)
            else
                0, 0
            end
            Δg::T = (imps_g_max - imps_g_min)
            if Δg > 0 
                g_imin::Int = searchsortednearest(ax, imps_g_min)
                g_imax::Int = searchsortednearest(ax, imps_g_max)
                n_grid_points::Int = g_imax - g_imin + 1
                n::Int = Int(round(Δg / stepsize))
                if n > 2 * n_grid_points
                    stepsize = Δg / (4 * n_grid_points)
                end
            else
                stepsize = 1
            end
            if iszero(stepsize) stepsize = 1 end
            push!(stepsizes, stepsize)
        end
        append!(samples, filter( x-> x in object.geometry, sample(g, stepsizes)) )
    end
    object_gridpoints = unique!([find_closest_gridpoint(sample_point, grid) for sample_point in samples])
    return object_gridpoints
end


function paint_object(det::SolidStateDetector{T}, object::AbstractObject{T}, grid::CylindricalGrid{T}, ::Val{:φ}, φ::T )  where {T <: SSDFloat}
    closest_φ_idx=searchsortednearest(grid.axes[2].ticks, φ)
    stepsize::Vector{T}= [minimum(diff(grid.axes[1].ticks)), IntervalSets.width(grid.axes[2].interval) == 0.0 ? 0.05236 : minimum(diff(grid.axes[2].ticks)), minimum(diff(grid.axes[3].ticks))]
    stepsize /= 2
    samples = filter(x-> x in object.geometry, vcat([sample(g, stepsize) for g in object.geometry_positive]...))
    object_gridpoints = unique!([find_closest_gridpoint(sample_point,grid) for sample_point in samples])
    return filter(x -> x[2]==closest_φ_idx, object_gridpoints)
end
function paint_object(det::SolidStateDetector{T}, object::AbstractObject{T}, grid::CylindricalGrid{T}, ::Val{:r}, r::T )  where {T <: SSDFloat}
    return CartesianPoint{T}[]
end
function paint_object(det::SolidStateDetector{T}, object::AbstractObject{T}, grid::Grid{T}, ::Val{:z}, z::T )  where {T <: SSDFloat}
    return CartesianPoint{T}[]
end
function paint_object(det::SolidStateDetector{T}, object::AbstractObject{T}, grid::CartesianGrid{T}, ::Val{:x}, x::T )  where {T <: SSDFloat}
    return CartesianPoint{T}[]
end
function paint_object(det::SolidStateDetector{T}, object::AbstractObject{T}, grid::CartesianGrid{T}, ::Val{:y}, y::T )  where {T <: SSDFloat}
    return CartesianPoint{T}[]
end
