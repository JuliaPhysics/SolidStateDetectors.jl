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
