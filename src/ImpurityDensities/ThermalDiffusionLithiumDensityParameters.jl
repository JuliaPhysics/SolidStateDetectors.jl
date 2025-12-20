function _ThermalDiffusionLithiumParameters(
    config::AbstractDict;
    T::Type = Float32,
    temperature_ranges::AbstractVector{<:AbstractDict} = config["annealing_temperature_ranges"] ,
    a::Union{Real, String} = config["experimental_parameters"]["a"],
    b::Union{Real, String} = config["experimental_parameters"]["b"],
    input_units::Union{Missing, NamedTuple} = missing
)

    temperature_unit = !ismissing(input_units) ? input_units.temperature : internal_temperature_unit
    diffusivity_unit = !ismissing(input_units) ?
        input_units.length^2 / internal_time_unit : internal_length_unit^2/internal_time_unit

temp_ranges_parameters = Vector{LithiumDiffusionParameters{T}}()

for r in temperature_ranges
    push!(temp_ranges_parameters,
        LithiumDiffusionParameters{T}(
            _parse_value(T, r["T_min"], temperature_unit),
            _parse_value(T, r["T_max"], temperature_unit),
            _parse_value(T, r["D0"], diffusivity_unit),
            _parse_value(T, r["H"], u"cal")
        )
    )
end

for i in 1:length(temp_ranges_parameters)
    if(temp_ranges_parameters[i].T_max < temp_ranges_parameters[i].T_min)
        throw(ConfigFileError("Invalid annealing temperature range $i: T_max must be ≥ T_min."))
    end
    if(i>1 && temp_ranges_parameters[i-1].T_max != temp_ranges_parameters[i].T_min)
        throw(ConfigFileError("Annealing temperature ranges must be contiguous and increasing: range $i-1 ends at a different temperature than range $i starts."))
    end
end

diffusion_params = Tuple(temp_ranges_parameters)

saturation_params = LithiumSaturationParameters{T}(
    T(a),
    T(b)
)

return ThermalDiffusionLithiumDensityParameters(
    diffusion_params,
    saturation_params
)
end

const minimum_num_of_temp_ranges = 2

function ThermalDiffusionLithiumParameters(config::AbstractDict, input_units::Union{Missing, NamedTuple} = missing; kwargs...)
    if !haskey(config, "annealing_temperature_ranges") 
        throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'annealing_temperature_ranges'.")) 
    elseif length(config["annealing_temperature_ranges"])<minimum_num_of_temp_ranges
        throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs at least $minimum_num_of_temp_ranges entries in 'annealing_temperature_ranges'."))
    end 

    for i in 1:length(config["annealing_temperature_ranges"])
        for axis in ("T_min", "T_max", "D0", "H")
            if !haskey(config["annealing_temperature_ranges"][i], axis)
                throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'annealing_temperature_ranges[$i]/$(axis)'."))
            end
        end
    end

    if !haskey(config, "experimental_parameters") 
        throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'experimental_parameters'.")) 
    elseif !haskey(config["experimental_parameters"], "a") 
        throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'experimental_parameters/a'.")) 
    elseif !haskey(config["experimental_parameters"], "b") 
        throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'experimental_parameters/b'.")) 
    end

    _ThermalDiffusionLithiumParameters(config; input_units, kwargs...)
end

const default_ThermalDiffusionLithiumDensity_config_file = joinpath(get_path_to_example_config_files(), "ThermalDiffusionLithiumDensity/Thermal_Diffusion_Config.yaml")
function ThermalDiffusionLithiumParameters(config_filename::AbstractString = default_ThermalDiffusionLithiumDensity_config_file, input_units::Union{Missing, NamedTuple} = missing; kwargs...)
    ThermalDiffusionLithiumParameters(parse_config_file(config_filename), input_units; kwargs...)
end