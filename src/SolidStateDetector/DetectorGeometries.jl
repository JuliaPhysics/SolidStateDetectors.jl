include("Object.jl")
include("Passive.jl")
include("Contacts.jl")
include("Semiconductor.jl")
include("AbstractVirtualVolume.jl")
include("TransitionLayer.jl")
include("SigGenInterface.jl")
include("SolidStateDetector.jl")


"""
    parse_config_file(filename::AbstractString)::SolidStateDetector{T} where {T <: SSDFloat}

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

function sample(c::SolidStateDetector{T}, ::Type{Cartesian}, sampling...)::Vector{CartesianPoint{T}} where {T <: SSDFloat}
    imp::Vector{CartesianPoint{T}} = vcat(
        CartesianPoint.(sample(c.semiconductor.geometry, sampling...)),
        [CartesianPoint.(sample(g.geometry, sampling...)) for object in skipmissing((c.contacts, c.passives)) for g in object]...)
    unique!(imp)
end

function sample(c::SolidStateDetector{T}, ::Type{Cylindrical}, sampling...)::Vector{CylindricalPoint{T}} where {T <: SSDFloat}
    imp::Vector{CylindricalPoint{T}} = vcat(
    CylindricalPoint.(sample(c.semiconductor.geometry, sampling...)),
    [CylindricalPoint.(sample(g.geometry, sampling...)) for object in skipmissing((c.contacts, c.passives)) for g in object]...)
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
        if in(searchsortednearest(grid, p), contact) || in(searchsortednearest(grid, p), contact) || p in contact || p in contact 
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

function get_ρ_and_ϵ(pt::AbstractCoordinatePoint{T}, ssd::SolidStateDetector{T}, medium::NamedTuple = material_properties[materials["vacuum"]])::Tuple{T, T, T} where {T <: SSDFloat}
    ρ_semiconductor::T = 0
    q_eff_fix::T = 0
    ϵ::T = medium.ϵ_r
    if pt in ssd.semiconductor
        ρ_semiconductor = get_charge_density(ssd.semiconductor, pt) 
        ϵ = ssd.semiconductor.material.ϵ_r
    elseif !ismissing(ssd.passives) && in(pt, ssd.passives)
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


function json_to_dict(inputfile::String)::Dict
    parsed_json_file = Dict()
    open(inputfile,"r") do f
        global parsed_json_file
        dicttext = readstring(f)
        parsed_json_file = JSON.parse(dicttext)
    end
    return parsed_json_file
end