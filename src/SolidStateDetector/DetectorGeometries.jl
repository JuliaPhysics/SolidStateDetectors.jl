include("Object.jl")
include("Passive.jl")
include("Contacts.jl")
include("Semiconductor.jl")
include("AbstractVirtualVolume.jl")
include("TransitionLayer.jl")
include("SigGenInterface.jl")
include("SolidStateDetector.jl")


"""
    SolidStateDetector{T}(filename::AbstractString)::SolidStateDetector{T} where {T <: SSDFloat}

Reads in a config file and returns an Detector struct which holds all information specified in the config file.
Currently supported formats for the config file: .json, .yaml
"""
function parse_config_file(filename::AbstractString)::Dict{Any,Any}
    if endswith(filename, ".toml")
        error("Currently only .json and .yaml files are supported.")
    elseif endswith(filename, ".json")
        dicttext = read(filename, String)
        parsed_dict = JSON.parse(dicttext)
        scan_and_merge_included_json_files!(parsed_dict, filename)
    elseif endswith(filename, ".yaml")
        parsed_dict = YAML.load_file(filename)
        scan_and_merge_included_json_files!(parsed_dict, filename)
    elseif endswith(filename, ".config")
        siggen_dict = readsiggen(filename)
        parsed_dict = siggentodict(siggen_dict)
    else
        error("Currently only .json and .yaml files are supported.")
    end
    parsed_dict
end

function yaml2json_convert(filename::String)
    data = YAML.load(open(filename), dicttype = DataStructures.OrderedDict)
    if isempty(data)
        @warn "The file $(filename) is empty and will not be converted."
    else
        open(replace(filename, ".yaml" => ".json"),"w") do f 
            JSON.print(f, data, 4) 
        end
    end
end
function yaml2json(directory::String)# or filename
    if endswith(directory, ".yaml")
        yaml2json_convert(directory)
    else
        yamlfiles = filter(x->endswith(x,".yaml"), joinpath.(directory, readdir(directory)))
        for filename in yamlfiles
            yaml2json_convert(filename)
        end
    end
end

function json2yaml_convert(filename::String)::Nothing
    data = JSON.parsefile(filename, dicttype = DataStructures.OrderedDict)
    if isempty(data)
        @warn "The file $(filename) is empty and will not be converted."
    else
        YAML.write_file(replace(filename, ".json" => ".yaml"), data)
    end

end
function json2yaml(directory::String)# or filename
    if endswith(directory, ".json")
        json2yaml_convert(directory)
    else
        jsonfiles = filter(x->endswith(x,".json"), joinpath.(directory, readdir(directory)))
        for filename in jsonfiles
            json2yaml_convert(filename)
        end
    end
end

function scan_and_merge_included_json_files!(parsed_dict, config_filename::AbstractString)
    key_word = "include"
    config_dir = dirname(config_filename)
    for k in keys(parsed_dict)
        is_subdict = typeof(parsed_dict[k]) <: Dict
        if !is_subdict && string(k) != key_word
            typeof(parsed_dict[k]) <: Array ? is_subdict = true : is_subdict = false
        end
        if is_subdict
            scan_and_merge_included_json_files!(parsed_dict[k], config_filename)
        elseif string(k) == key_word
            files = []
            if typeof(parsed_dict[k]) <: Array
                append!(files, parsed_dict[k])
            else
                push!(files, parsed_dict[k])
            end
            for file in files
                filepath = joinpath(config_dir, file)
                if isfile(filepath)
                    tmp = parse_config_file(filepath)
                    for sub_k in keys(tmp)
                        parsed_dict[sub_k] = tmp[sub_k]
                    end
                end
            end
            delete!(parsed_dict, k)
        end
    end
end


function SolidStateDetector{T}(filename::AbstractString)::SolidStateDetector{T} where {T <: SSDFloat}
    parsed_dict = parse_config_file(filename)
    return SolidStateDetector{T}(parsed_dict)
end

# function SolidStateDetector(T::Type{<:AbstractFloat} = Float32, filename::AbstractString = SSD_examples[:InvertedCoax])::SolidStateDetector{T}
#     SolidStateDetector{T}(filename)
# end

function SolidStateDetector(filename::AbstractString)::SolidStateDetector{Float32}
    SolidStateDetector{Float32}(filename)
end

function sample(c::SolidStateDetector{T, Cartesian}, sampling...)::Vector{CartesianPoint{T}} where {T <: SSDFloat}
    imp::Vector{CartesianPoint{T}} = vcat(
        CartesianPoint.(sample(c.semiconductor.geometry, sampling...)),
        [CartesianPoint.(sample(g.geometry, sampling...)) for object in (c.contacts, c.passives) for g in object]...)
    unique!(imp)
end

function sample(c::SolidStateDetector{T, Cylindrical}, sampling...)::Vector{CylindricalPoint{T}} where {T <: SSDFloat}
    imp::Vector{CylindricalPoint{T}} = vcat(
    CylindricalPoint.(sample(c.semiconductor.geometry, sampling...)),
    [CylindricalPoint.(sample(g.geometry, sampling...)) for object in (c.contacts, c.passives) for g in object]...)
    unique!(imp)
end


function is_boundary_point(c::SolidStateDetector, pt::AbstractCoordinatePoint{T}, ax1::Vector{T}, ax2::Vector{T}, ax3::Vector{T})::Tuple{Bool, Real, Int} where {T <: SSDFloat}
    # if false #!(p in c)
    #     return false, T(0), 0
    # else
    for contact in c.contacts
        if in(pt, contact, ax1)
            return true, contact.potential, contact.id
        end
    end
    for external_part in c.external_parts
        if in(pt, external_part, ax1)
            return true, external_part.potential, external_part.id
        end
    end
    return false, 0.0, 0
    # end
end


function point_type(c::SolidStateDetector{T}, grid::Grid{T, 3}, p::CylindricalPoint{T})::Tuple{UInt8, Int, CartesianVector{T}} where {T <: SSDFloat}
    surface_normal::CartesianVector{T} = CartesianVector{T}(0, 0, 0) # need undef version for this
    for contact in c.contacts
        if in(searchsortednearest(grid, geom_round(p)), contact) || in(searchsortednearest(grid, p), contact) || geom_round(p) in contact || p in contact 
            return CD_ELECTRODE::UInt8, contact.id, surface_normal
        end
    end
    on_surface, surface_normal = is_surface_point_and_normal_vector(c, p)
    if on_surface
        for contact in c.contacts
            if in(searchsortednearest(grid, p), contact) #&& abs(sum(sp[2])) > 1
                return CD_ELECTRODE::UInt8, contact.id, surface_normal
            else
                return CD_FLOATING_BOUNDARY::UInt8, -1, surface_normal
            end
        end
    elseif !(p in c)
        return CD_OUTSIDE::UInt8, -1, surface_normal
    else
        return CD_BULK::UInt8, -1, surface_normal
    end
end


# Point types for charge drift
const CD_ELECTRODE = 0x00
const CD_OUTSIDE = 0x01
const CD_BULK = 0x02
const CD_FLOATING_BOUNDARY = 0x04 # not 0x03, so that one could use bit operations here...

"""
    For charge drift...
"""
function point_type(c::SolidStateDetector{T}, grid::Grid{T, 3}, p::CartesianPoint{T})::Tuple{UInt8, Int, CartesianVector{T}} where {T <: SSDFloat}
    surface_normal::CartesianVector{T} = CartesianVector{T}(0, 0, 0) # need undef version for this
    for contact in c.contacts
        if p in contact #|| geom_round(p) in contact #|| in(go_to_nearest_gridpoint(c,p), contact, c.rs) || in(go_to_nearest_gridpoint(c,geom_round(p)), contact, c.rs)
            return CD_ELECTRODE::UInt8, contact.id, surface_normal
        end
    end
    on_surface, surface_normal = is_surface_point_and_normal_vector(c, p) # surface_normal::CartesianVector{T}
    if on_surface
        for contact in c.contacts
            if in(searchsortednearest(grid, p), contact)
                return CD_ELECTRODE::UInt8, contact.id, surface_normal
            else
                return CD_FLOATING_BOUNDARY::UInt8, -1, surface_normal
            end
        end
    elseif !(p in c)
        return CD_OUTSIDE::UInt8, -1, surface_normal
    else
        return CD_BULK::UInt8, -1, surface_normal
    end
end

# go_to_nearest_gridpoint
function searchsortednearest(grid::Grid{T, 3, Cylindrical}, pt::CylindricalPoint{T})::CylindricalPoint{T} where {T <: SSDFloat}
    idx1::Int = searchsortednearest(grid.axes[1].ticks, pt.r)
    idx2::Int = searchsortednearest(grid.axes[2].ticks, pt.φ)
    idx3::Int = searchsortednearest(grid.axes[3].ticks, pt.z)
    CylindricalPoint{T}(grid.axes[1].ticks[idx1], grid.axes[2].ticks[idx2], grid.axes[3].ticks[idx3])
end
function searchsortednearest(grid::Grid{T, 3, Cartesian}, pt::CartesianPoint{T})::CartesianPoint{T} where {T <: SSDFloat}
    idx1::Int = searchsortednearest(grid.axes[1].ticks, pt.x)
    idx2::Int = searchsortednearest(grid.axes[2].ticks, pt.y)
    idx3::Int = searchsortednearest(grid.axes[3].ticks, pt.z)
    CartesianPoint{T}(grid.axes[1].ticks[idx1], grid.axes[2].ticks[idx2], grid.axes[3].ticks[idx3])
end

function is_surface_point_and_normal_vector(c::SolidStateDetector{T}, p::CylindricalPoint{T})::Tuple{Bool, CartesianVector{T}} where {T <: SSDFloat}
    if !(p in c) # contacts are already checked in 
        return false, CartesianPoint{T}(0, 0, 0)
    end
    n::MVector{3,T} = @MVector T[0, 0, 0]
    look_around::Vector{Bool} = [   CylindricalPoint{T}(prevfloat(p.r), p.φ, p.z) in c,
                                    CylindricalPoint{T}(nextfloat(p.r), p.φ, p.z) in c,
                                    CylindricalPoint{T}(p.r, prevfloat(p.φ), p.z) in c,
                                    CylindricalPoint{T}(p.r, nextfloat(p.φ), p.z) in c,
                                    CylindricalPoint{T}(p.r, p.φ, prevfloat(p.z)) in c,
                                    CylindricalPoint{T}(p.r, p.φ, nextfloat(p.z)) in c]
    if !(false in look_around)
        return false , CartesianPoint{T}(n...)
    else
        if (look_around[1]==false) n[1] -= 1 end
        if (look_around[2]==false) n[1] += 1 end
        if (look_around[3]==false) n[2] -= 1 end
        if (look_around[4]==false) n[2] += 1 end
        if (look_around[5]==false) n[3] -= 1 end
        if (look_around[6]==false) n[3] += 1 end
        # println(look_around , " " , n)
        Rα::SMatrix{3,3,T} = @SArray([cos(p.φ) -1*sin(p.φ) 0;sin(p.φ) cos(p.φ) 0;0 0 1])
        # return true, geom_round(CartesianVector{T}((Rα * n)...))
        return true, CartesianVector{T}((Rα * n)...)
    end
end
function is_surface_point_and_normal_vector(c::SolidStateDetector{T}, p::CartesianPoint{T})::Tuple{Bool, CartesianVector{T}} where T
    if !(p in c) 
        return false, CartesianPoint{T}(0, 0, 0)
    end
    n::MVector{3,T} = @MVector T[0.0,0.0,0.0]
    look_around::Vector{Bool} = [   CartesianPoint{T}(prevfloat(p.x), p.y, p.z) in c,
                                    CartesianPoint{T}(nextfloat(p.x), p.y, p.z) in c,
                                    CartesianPoint{T}(p.x, prevfloat(p.y), p.z) in c,
                                    CartesianPoint{T}(p.x, nextfloat(p.y), p.z) in c,
                                    CartesianPoint{T}(p.x, p.y, prevfloat(p.z)) in c,
                                    CartesianPoint{T}(p.x, p.y, nextfloat(p.z)) in c]
    if !(false in look_around)
        return false, n
    else
        if (look_around[1] == false) n[1] -= 1 end
        if (look_around[2] == false) n[1] += 1 end
        if (look_around[3] == false) n[2] -= 1 end
        if (look_around[4] == false) n[2] += 1 end
        if (look_around[5] == false) n[3] -= 1 end
        if (look_around[6] == false) n[3] += 1 end
        # println(look_around , " " , n)
        return true, n
    end
end


function get_charge_density(sc::Semiconductor{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    get_impurity_density(sc.impurity_density_model, pt) * elementary_charge
end
function get_charge_density(p::Passive{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    get_charge_density(p.charge_density_model, pt)
end

function get_ρ_and_ϵ(pt::AbstractCoordinatePoint{T}, ssd::SolidStateDetector{T})::Tuple{T, T, T} where {T <: SSDFloat}
    ρ_semiconductor::T = 0
    q_eff_fix::T = 0
    ϵ::T = ssd.medium.ϵ_r
    if pt in ssd.semiconductor
        ρ_semiconductor = get_charge_density(ssd.semiconductor, pt) 
        ϵ = ssd.semiconductor.material.ϵ_r
    elseif in(pt, ssd.passives)
        for ep in ssd.passives
            if pt in ep
                q_eff_fix = get_charge_density(ep, pt)
                ϵ = ep.material.ϵ_r
                break
            end
        end
    end
    return ρ_semiconductor, ϵ, q_eff_fix
end


function set_pointtypes_and_fixed_potentials!(pointtypes::Array{PointType, N}, potential::Array{T, N},
        grid::Grid{T, N, Cylindrical}, ssd::SolidStateDetector{T}; weighting_potential_contact_id::Union{Missing, Int} = missing)::Nothing where {T <: SSDFloat, N}

    channels::Array{Int, 1} = if !ismissing(weighting_potential_contact_id)
        [weighting_potential_contact_id]
    else
        Int[]
    end

    axr::Vector{T} = grid.axes[1].ticks
    axφ::Vector{T} = grid.axes[2].ticks
    axz::Vector{T} = grid.axes[3].ticks

    for iz in axes(potential, 3)
        z::T = axz[iz]
        for iφ in axes(potential, 2)
            φ::T = axφ[iφ]
            for ir in axes(potential, 1)
                r::T = axr[ir]
                pt::CylindricalPoint{T} = CylindricalPoint{T}( r, φ, z )

                for passive in ssd.passives
                    if passive.potential != :floating
                        if pt in passive
                            potential[ ir, iφ, iz ] = if ismissing(weighting_potential_contact_id)
                                passive.potential
                            else
                                0
                            end
                            pointtypes[ ir, iφ, iz ] = zero(PointType)
                        end
                    end
                end
                if in(pt, ssd)
                    pointtypes[ ir, iφ, iz ] += pn_junction_bit
                end
                for contact in ssd.contacts
                    if pt in contact
                        potential[ ir, iφ, iz ] = if ismissing(weighting_potential_contact_id)
                            contact.potential
                        else
                            contact.id == weighting_potential_contact_id ? 1 : 0
                        end
                        pointtypes[ ir, iφ, iz ] = zero(PointType)
                    end
                end
            end
        end
    end

    for contact in ssd.contacts
        pot::T = if ismissing(weighting_potential_contact_id)
            contact.potential
        else
            contact.id == weighting_potential_contact_id ? 1 : 0
        end
        contact_gridpoints = paint_object(contact, grid)
        for gridpoint in contact_gridpoints
            potential[ gridpoint... ] = pot
            pointtypes[ gridpoint... ] = zero(PointType)
        end
    end
    nothing
end

function set_pointtypes_and_fixed_potentials!(pointtypes::Array{PointType, N}, potential::Array{T, N},
    grid::Grid{T, N, Cartesian}, ssd::SolidStateDetector{T}; weighting_potential_contact_id::Union{Missing, Int} = missing)::Nothing where {T <: SSDFloat, N}

    channels::Array{Int, 1} = if !ismissing(weighting_potential_contact_id)
        [weighting_potential_contact_id]
    else
        Int[]
    end

    axx::Vector{T} = grid.axes[1].ticks
    axy::Vector{T} = grid.axes[2].ticks
    axz::Vector{T} = grid.axes[3].ticks

    for iz in axes(potential, 3)
        z::T = axz[iz]
        for iy in axes(potential, 2)
            y::T = axy[iy]
            for ix in axes(potential, 1)
                x::T = axx[ix]
                pt::CartesianPoint{T} = CartesianPoint{T}( x, y, z )

                for passive in ssd.passives
                    if passive.potential != :floating
                        if pt in passive
                            potential[ ix, iy, iz ] = if ismissing(weighting_potential_contact_id)
                                passive.potential
                            else
                                0
                            end
                            pointtypes[ ix, iy, iz ] = zero(PointType)
                        end
                    end
                end
                if in(pt, ssd)
                    pointtypes[ ix, iy, iz ] += pn_junction_bit
                end
                for contact in ssd.contacts
                    if pt in contact
                        potential[ ix, iy, iz ] = if ismissing(weighting_potential_contact_id)
                            contact.potential
                        else
                            contact.id == weighting_potential_contact_id ? 1 : 0
                        end
                        pointtypes[ ix, iy, iz ] = zero(PointType)
                    end
                end
            end
        end
    end
    for contact in ssd.contacts
        pot::T = if ismissing(weighting_potential_contact_id)
            contact.potential
        else
            contact.id == weighting_potential_contact_id ? 1 : 0
        end
        contact_gridpoints = paint_object(contact, grid)
        for gridpoint in contact_gridpoints
            potential[ gridpoint... ] = pot
            pointtypes[ gridpoint... ] = zero(PointType)
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

function Grid(  detector::SolidStateDetector{T, Cylindrical};
                init_grid_size::Union{Missing, NTuple{3, Int}} = missing,
                init_grid_spacing::Union{Missing, Tuple{<:Real,<:Real,<:Real}} = missing,
                for_weighting_potential::Bool = false,
                min_n_ticks::Int = 10,
                full_2π::Bool = false)::CylindricalGrid{T} where {T}
    if ismissing(init_grid_size)
        world_diffs = [(getproperty.(detector.world.intervals, :right) .- getproperty.(detector.world.intervals, :left))...]
        world_diffs[2] = world_diffs[2] * 0.3 * detector.world.intervals[1].right # in radiance
        inds::Vector{Int} = sortperm([world_diffs...])
        ratio::T = min_n_ticks * if world_diffs[inds[1]] > 0
            inv(world_diffs[inds[1]])
        elseif world_diffs[inds[2]] > 0
            inv(world_diffs[inds[2]])
        elseif world_diffs[inds[3]] > 0
            inv(world_diffs[inds[3]])
        else
            error("This should not happen... World has no dimension")
        end
        init_grid_size_1::Int = convert(Int, round(ratio * world_diffs[1], RoundUp))
        init_grid_size_2::Int = convert(Int, round(ratio * world_diffs[2], RoundUp))
        init_grid_size_3::Int = convert(Int, round(ratio * world_diffs[3], RoundUp))
        init_grid_size::NTuple{3, Int} = NTuple{3, T}( [init_grid_size_1, init_grid_size_2, init_grid_size_3] )
    end

    init_grid_spacing, use_spacing::Bool = if !ismissing(init_grid_spacing)
        T.(init_grid_spacing), true
    else
        missing, false
    end
    
    samples::Vector{CylindricalPoint{T}} = sample(detector)
   
    important_r_points::Vector{T} = map(p -> p.r, samples)
    important_φ_points::Vector{T} = map(p -> p.φ, samples)
    important_z_points::Vector{T} = map(p -> p.z, samples)

    push!(important_r_points, detector.world.intervals[1].left)
    push!(important_r_points, detector.world.intervals[1].right)
    important_r_points = unique!(sort!(geom_round.(important_r_points)))
    push!(important_z_points, detector.world.intervals[3].left)
    push!(important_z_points, detector.world.intervals[3].right)
    important_z_points = unique!(sort!(geom_round.(important_z_points)))
    push!(important_φ_points, detector.world.intervals[2].left)
    push!(important_φ_points, detector.world.intervals[2].right)
    important_φ_points = unique!(sort!(geom_round.(important_φ_points)))

    # r
    L, R, BL, BR = get_boundary_types(detector.world.intervals[1])
    int_r = Interval{L, R, T}(detector.world.intervals[1].left, detector.world.intervals[1].right)
    ax_r::DiscreteAxis{T, BL, BR} = if use_spacing
        DiscreteAxis{BL, BR}(int_r, step = init_grid_spacing[1])
    else
        DiscreteAxis{BL, BR}(int_r, length = init_grid_size[1])
    end
    rticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_r, important_r_points, atol = minimum(diff(ax_r.ticks))/4)
    ax_r = DiscreteAxis{T, BL, BR}(int_r, rticks)


    # φ
    L, R, BL, BR = get_boundary_types(detector.world.intervals[2])
    int_φ = Interval{L, R, T}(detector.world.intervals[2].left, detector.world.intervals[2].right)
    if full_2π == true || (for_weighting_potential && (detector.world.intervals[2].left != detector.world.intervals[2].right))
        L, R, BL, BR = :closed, :open, :periodic, :periodic
        int_φ = Interval{L, R, T}(0, 2π)
    end
    ax_φ = if int_φ.left == int_φ.right
        DiscreteAxis{T, BL, BR}(int_φ, T[int_φ.left])
    else
        if use_spacing
            DiscreteAxis{BL, BR}(int_φ, step = init_grid_spacing[2])
        else
            DiscreteAxis{BL, BR}(int_φ, length = init_grid_size[2])
        end
    end
    if length(ax_φ) > 1
        φticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_φ, important_φ_points, atol = minimum(diff(ax_φ.ticks))/4)
        ax_φ = typeof(ax_φ)(int_φ, φticks)
    end
    if isodd(length(ax_φ)) && length(ax_φ) > 1 # must be even
        int_φ = ax_φ.interval
        φticks = ax_φ.ticks
        push!(φticks, geom_round((φticks[end] + φticks[end-1]) * 0.5))
        sort!(φticks)
        ax_φ = typeof(ax_φ)(int_φ, φticks) # must be even
    end
    if length(ax_φ) > 1
        @assert iseven(length(ax_φ)) "CylindricalGrid must have even number of points in φ."
    end

    #z
    L, R, BL, BR = get_boundary_types(detector.world.intervals[3])
    int_z = Interval{L, R, T}(detector.world.intervals[3].left, detector.world.intervals[3].right)
    ax_z::DiscreteAxis{T, BL, BR} = if use_spacing
        DiscreteAxis{BL, BR}(int_z, step = init_grid_spacing[3])
    else
        DiscreteAxis{BL, BR}(int_z, length = init_grid_size[3])
    end
    zticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_z, important_z_points, atol=minimum(diff(ax_z.ticks))/2)
    ax_z = typeof(ax_z)(int_z, zticks)
    if isodd(length(ax_z)) # must be even
        int_z = ax_z.interval
        zticks = ax_z.ticks
        push!(zticks, geom_round((zticks[end] + zticks[end-1]) * 0.5))
        sort!(zticks)
        ax_z = typeof(ax_z)(int_z, zticks) # must be even
    end
    @assert iseven(length(ax_z)) "CylindricalGrid must have even number of points in z."

    return CylindricalGrid{T}( (ax_r, ax_φ, ax_z) )
end


function Grid(  detector::SolidStateDetector{T, Cartesian};
                init_grid_size::Union{Missing, NTuple{3, Int}} = missing,
                init_grid_spacing::Union{Missing, Tuple{<:Real,<:Real,<:Real,}} = missing,
                min_n_ticks::Int = 10,
                for_weighting_potential::Bool = false)::CartesianGrid3D{T} where {T}

    if ismissing(init_grid_size)
        world_diffs = [(getproperty.(detector.world.intervals, :right) .- getproperty.(detector.world.intervals, :left))...]
        inds::Vector{Int} = sortperm([world_diffs...])
        ratio::T = min_n_ticks * if world_diffs[inds[1]] > 0
            inv(world_diffs[inds[1]])
        elseif world_diffs[inds[2]] > 0
            inv(world_diffs[inds[2]])
        elseif world_diffs[inds[3]] > 0
            inv(world_diffs[inds[3]])
        else
            error("This should not happen... World has no dimension")
        end
        init_grid_size_1::Int = convert(Int, round(ratio * world_diffs[1], RoundUp))
        init_grid_size_2::Int = convert(Int, round(ratio * world_diffs[2], RoundUp))
        init_grid_size_3::Int = convert(Int, round(ratio * world_diffs[3], RoundUp))
        init_grid_size::NTuple{3, Int} = NTuple{3, T}( [init_grid_size_1, init_grid_size_2, init_grid_size_3] )
    end

    samples::Vector{CartesianPoint{T}} = sample(detector)
    
    important_x_points::Vector{T} = geom_round.(map(p -> p.x, samples))
    important_y_points::Vector{T} = geom_round.(map(p -> p.y, samples))
    important_z_points::Vector{T} = geom_round.(map(p -> p.z, samples))

    init_grid_spacing, use_spacing::Bool = if !ismissing(init_grid_spacing)
        T.(init_grid_spacing), true
    else
        missing, false
    end

    # x
    L, R, BL, BR = get_boundary_types(detector.world.intervals[1])
    int_x = Interval{L, R, T}(detector.world.intervals[1].left, detector.world.intervals[1].right)
    ax_x::DiscreteAxis{T, BL, BR} = if use_spacing
        DiscreteAxis{BL, BR}(int_x, step = init_grid_spacing[1])
    else
        DiscreteAxis{BL, BR}(int_x, length = init_grid_size[1])
    end
    xticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_x, important_x_points, atol = minimum(diff(ax_x.ticks)) / 2)
    ax_x = typeof(ax_x)(int_x, xticks)
    if isodd(length(ax_x)) # RedBlack dimension must be of even length
        xticks = ax_x.ticks
        push!(xticks, geom_round((xticks[end] + xticks[end-1]) * 0.5))
        sort!(xticks)
        ax_x = DiscreteAxis{T, BL, BR}(int_x, xticks) # must be even
    end
    @assert iseven(length(ax_x)) "CartesianGrid3D must have even number of points in z."

    # y
    L, R, BL, BR = get_boundary_types(detector.world.intervals[2])
    int_y = Interval{L, R, T}(detector.world.intervals[2].left, detector.world.intervals[2].right)
    ax_y::DiscreteAxis{T, BL, BR} = if use_spacing
        DiscreteAxis{BL, BR}(int_y, step = init_grid_spacing[2])
    else
        DiscreteAxis{BL, BR}(int_y, length = init_grid_size[2])
    end
    yticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_y, important_y_points, atol = minimum(diff(ax_y.ticks)) / 2)
    ax_y = typeof(ax_y)(int_y, yticks)

    # z
    L, R, BL, BR = get_boundary_types(detector.world.intervals[3])
    int_z = Interval{L, R, T}(detector.world.intervals[3].left, detector.world.intervals[3].right)
    ax_z::DiscreteAxis{T, BL, BR} = if use_spacing
        DiscreteAxis{BL, BR}(int_z, step = init_grid_spacing[3])
    else
        DiscreteAxis{BL, BR}(int_z, length = init_grid_size[3])
    end
    zticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_z, important_z_points, atol = minimum(diff(ax_z.ticks)) / 2)
    ax_z = typeof(ax_z)(int_z, zticks)

    return CartesianGrid3D{T}( (ax_x, ax_y, ax_z) )
end