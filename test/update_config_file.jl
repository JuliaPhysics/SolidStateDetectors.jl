using SolidStateDetectors.ConstructiveSolidGeometry: CSGUnion, CSGIntersection, CSGDifference,
        internal_unit_length, internal_unit_angle, CSG_dict, geom_round,
        parse_r_of_primitive, parse_height_of_primitive, parse_translate_vector, parse_rotation_matrix
using SolidStateDetectors: construct_units
import SolidStateDetectors.ConstructiveSolidGeometry: Geometry, Dictionary
using SolidStateDetectors.ConstructiveSolidGeometry: AbstractGeometry, TranslatedGeometry, ScaledGeometry, RotatedGeometry, 
    CSGUnion, CSGDifference, CSGIntersection, Cone, CartesianVector, HexagonalPrism
using DataStructures: OrderedDict
using IntervalSets
using Rotations
using Unitful
using JSON
using YAML


#old implementation of Geometry:    
function Geometry(::Type{T}, t::Type{CSGUnion}, dict::Union{Dict, OrderedDict}, input_units::NamedTuple) where {T}
    @assert haskey(dict, "parts") "Please specify 'parts' of the '$(dict["type"])'."
    sum( map(x-> Geometry(T, x, input_units), dict["parts"]) )
end

function Geometry(::Type{T}, ::Type{CSGIntersection}, dict::Union{Dict, OrderedDict}, input_units::NamedTuple) where {T}
    @assert haskey(dict, "parts") "Please specify 'parts' of the '$(dict["type"])'."
    parts = map(x-> Geometry(T, x, input_units), dict["parts"]) 
    reduce(&, parts)
end

function Geometry(::Type{T}, ::Type{CSGDifference}, dict::Union{Dict, OrderedDict}, input_units::NamedTuple) where {T}
    @assert haskey(dict, "parts") "Please specify 'parts' of the '$(dict["type"])'."
    Geometry(T, dict["parts"][1], input_units) - sum( map(x-> Geometry(T, x, input_units), dict["parts"][2:end]) )
end

function Geometry(::Type{T}, dict::Union{Dict, OrderedDict}, input_units::NamedTuple) where {T}
    if haskey(dict, "translate")
        length_unit = input_units.length
        t::CartesianVector{T} = parse_translate_vector(T, dict["translate"], length_unit)
        gdict::Dict{Any,Any} = filter(p -> first(p) != "translate", dict)
        return Geometry(T, gdict, input_units) + t
    end
    Geometry(T, CSG_dict[dict["type"]], dict, input_units)
end

function Geometry(::Type{T}, filename::String, input_units::NamedTuple) where {T}
    @assert isfile(filename) "The given filename '$(filename)' does not lead to a valid file."
    dict = if endswith(filename, ".json")
        JSON.parsefile(filename)
    elseif endswith(filename, ".yaml")
        YAML.load_file(filename)
    else
        @error "Only JSON and YAML formats are supported at the moment."
    end
    Geometry(T, dict, input_units)
end


# UPDATE DUE TO BREAKING CHANGES TO TUBE/CONE: CENTERED AROUND ZERO
# AND REDEFINITION OF HEXAGONALPRISM
function update_primitives!(dict::AbstractDict, input_units::NamedTuple)
    dict_keys = keys(dict)
    #Translate all Tubes and Cone that are defined via "h" and where h is not zero
    if "type" in dict_keys && dict["type"] == "HexagonalPrism"
        rotZ = haskey(dict, "rotZ") ? dict["rotZ"] : 0
        rotZ -= input_units.angle == u"°" ? 30 : π/6
        rotZ = ustrip(uconvert(u"rad", rotZ * input_units.angle))
        delete!(dict, "rotZ")
        if rotZ != 0 dict["rotate"] = Dict("Z" => rotZ) end
    elseif "type" in dict_keys && dict["type"] in ["tube", "cone"] && "h" in dict_keys && dict["h"] != 0
        #add translate vector with z=0 if it does not exist
        if !("translate" in dict_keys) dict["translate"] = Dict{String,Any}("z" => 0) end
        dict["translate"]["z"] = geom_round(dict["translate"]["z"] + Float64(dict["h"] / 2))
        if dict["translate"]["z"] == 0 delete!(dict["translate"], "z") end
        if length(keys(dict["translate"])) == 0 delete!(dict, "translate") end
    end
    for k in dict_keys
        if typeof(dict[k]) <: AbstractDict
            dict[k] = update_primitives!(dict[k], input_units)
        elseif typeof(dict[k]) <: Array
            for i in eachindex(dict[k])
                if typeof(dict[k][i]) <: AbstractDict
                    dict[k][i] = update_primitives!(dict[k][i], input_units)
                end
            end
        end
    end
    dict
end

function update_geometry!(g::AbstractDict, T::DataType, input_units = (length = NoUnits, angle = NoUnits))
    g["geometry"] = Dictionary(Geometry(T, g["geometry"], input_units))
    return g
end

# HexagonalPrism (#146)
function Geometry(::Type{T}, ::Type{HexagonalPrism}, dict::AbstractDict, input_units::NamedTuple) where {T}
    length_unit = input_units.length
    r = parse_r_of_primitive(T, dict, length_unit)
    z = parse_height_of_primitive(T, dict, length_unit)
    rot = haskey(dict, "rotate") ? parse_rotation_matrix(T, dict["rotate"], u"rad") : one(RotMatrix{3,T})
    return RotatedGeometry(HexagonalPrism(T, r, z), rot)
end

function Dictionary(g::RotatedGeometry{T}) where {T}
    dict = Dictionary(g.p)
    mat = RotXYZ(inv(g.inv_r))
    if mat.theta1 == 0 && mat.theta2 == 0
        if mat.theta3 != 0 
            dict["Z"] = geom_round(HexagonalPrismAngleUnit == u"°" ? rad2deg(mat.theta3) : mat.theta3)
        end
    else
        dict["M"] = inv(g.inv_r)[:]
    end
    OrderedDict{String,Any}("rotate" => dict)
end


function restructure_config_file_dict!(config_file_dict::AbstractDict, T::DataType = Float64)

    input_units = construct_units(config_file_dict)
    update_primitives!(config_file_dict, input_units)
    global HexagonalPrismAngleUnit = input_units.angle

    t = broadcast(obj -> obj["type"], config_file_dict["objects"])
    
    semiconductors = config_file_dict["objects"][findall(t .== "semiconductor")]
    @assert length(semiconductors) == 1 "Config files need exactly one semiconductor per detector. This one has $(length(semiconductors)) semiconductors."
    semiconductor = semiconductors[1]
    
    semiconductor = begin 
        delete!(semiconductor, "type")
        if "charge_density_model" in keys(semiconductor)
            semiconductor["impurity_density"] = semiconductor["charge_density_model"]
            delete!(semiconductor, "charge_density_model")
        end
        update_geometry!(semiconductor, T)
        semiconductor
    end

    contacts = [
        begin 
        delete!(contact, "type")
        update_geometry!(contact, T)
        contact
        end
        for contact in config_file_dict["objects"][findall(t .== "contact")]
    ]

    passives = [
        begin 
        delete!(passive, "type")
        update_geometry!(passive, T)
        passive
        end
        for passive in config_file_dict["objects"][findall(t .== "passive")]
    ]

    virtual_drift_volumes = [
        begin 
        delete!(virtual_volume, "type")
        update_geometry!(virtual_volume, T)
        virtual_volume
        end
        for virtual_volume in config_file_dict["objects"][findall(t .== "virtual_drift_volume")]
    ]

    config_file_dict["detectors"] = []
    push!(config_file_dict["detectors"], Dict("bulk" => semiconductor) )
    if length(contacts) > 0 push!(config_file_dict["detectors"], Dict("contacts" => contacts)) end
    if length(passives) > 0 push!(config_file_dict["detectors"], Dict("passives" => passives)) end
    if length(virtual_drift_volumes) > 0 push!(config_file_dict["detectors"], Dict("virtual_drift_volumes" => virtual_drift_volumes)) end
    
    config_file_dict
    
end
            

"""
    update_config_file(filename::String)::String

Updates the format of config files from v0.5 to v0.6.
Currently supported formats for the config file: .json, .yaml.
Returns the new filename.
"""
function update_config_file(filename::String)::String
    @assert isfile(filename) "'$(filename)' does not exist!"
    
    if endswith(filename, ".json")
        data = JSON.parsefile(filename, dicttype = OrderedDict)
        restructure_config_file_dict!(data)
        new_filename = replace(filename, ".json" => "_updated.json")
        open(new_filename,"w") do f 
            JSON.print(f, data, 4) 
        end
    elseif endswith(filename, ".yaml")
        data = YAML.load_file(filename, dicttype = OrderedDict)
        restructure_config_file_dict!(data)
        new_filename = replace(filename, ".yaml" => "_updated.yaml")
        YAML.write_file(new_filename, data)
    else 
        error("Currently only .json and .yaml files are supported.")
        new_filename = ""
    end
    new_filename
end

