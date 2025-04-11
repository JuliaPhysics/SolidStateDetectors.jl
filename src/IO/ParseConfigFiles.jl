"""
    parse_config_file(filename::AbstractString; dicttype::Type = Dict)::dicttype

Reads in a configuration file and returns a parsed dictionary which holds all the
information specified in the configuration file.

Find detailed information on configuration files in [Configuration Files](@ref).

## Arguments 
* `filename::AbstractString`: File name of the configuration file. If the file is not
    in the same directory, a path to the file is required.

## Keywords
* `dicttype::Type`: Type of the dictionary returned.

!!! note
    Currently supported formats for the configuration files:
        - [YAML](https://github.com/JuliaData/YAML.jl): `filename` ends with `.yaml`.
        - [JSON](https://github.com/JuliaIO/JSON.jl): `filename` ends with `.json`.
        - [SigGen](https://github.com/radforddc/icpc_siggen): `filename` ends with `.config`.
"""
function parse_config_file(filename::AbstractString; dicttype::Type = Dict)::dicttype
    if endswith(filename, ".json")
        dicttext = read(filename, String)
        dict = JSON.parse(dicttext, dicttype = dicttype)
        scan_and_merge_included_json_files!(dict, filename)
    elseif endswith(filename, ".yaml")
        dict = YAML.load_file(filename, dicttype = dicttype)
        scan_and_merge_included_json_files!(dict, filename)
    elseif endswith(filename, ".config")
        siggen_dict = readsiggen(filename, dicttype = dicttype)
        dict = siggentodict(siggen_dict, dicttype = dicttype)
    else
        error("Currently only .json and .yaml files are supported.")
    end
    dict
end

function yaml2json_convert(filename::String)
    data = YAML.load(open(filename), dicttype = OrderedCollections.OrderedDict)
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
    data = JSON.parsefile(filename, dicttype = OrderedCollections.OrderedDict)
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

function scan_and_merge_included_json_files!(dict, config_filename::AbstractString)
    key_word = "include"
    config_dir = dirname(config_filename)
    for k in keys(dict)
        is_subdict = typeof(dict[k]) <: AbstractDict
        if !is_subdict && string(k) != key_word
            is_subdict = (typeof(dict[k]) <: AbstractArray)
        end
        if is_subdict
            scan_and_merge_included_json_files!(dict[k], config_filename)
        elseif string(k) == key_word
            files = []
            if typeof(dict[k]) <: AbstractArray
                append!(files, dict[k])
            else
                push!(files, dict[k])
            end
            for file in files
                filepath = joinpath(config_dir, file)
                if isfile(filepath)
                    tmp = parse_config_file(filepath)
                    for sub_k in keys(tmp)
                        dict[sub_k] = tmp[sub_k]
                    end
                end
            end
            delete!(dict, k)
        end
    end
end
