


function _ThermalDiffusionLithiumParameters(
    config::AbstractDict;
    T::Type = Float32,
    lowT_Tmin::Union{RealQuantity, String} =
        config["impurity_density_profile"]["annealing_temperature_ranges"]["low_T"]["T_min"],
    lowT_Tmax::Union{RealQuantity, String} =
        config["impurity_density_profile"]["annealing_temperature_ranges"]["low_T"]["T_max"],
    lowT_D0::Union{RealQuantity, String} =
        config["impurity_density_profile"]["annealing_temperature_ranges"]["low_T"]["D0"],
    lowT_H::Union{RealQuantity, String} =
        config["impurity_density_profile"]["annealing_temperature_ranges"]["low_T"]["H"],
    highT_Tmin::Union{RealQuantity, String} =
        config["impurity_density_profile"]["annealing_temperature_ranges"]["high_T"]["T_min"],
    highT_Tmax::Union{RealQuantity, String} =
        config["impurity_density_profile"]["annealing_temperature_ranges"]["high_T"]["T_max"],
    highT_D0::Union{RealQuantity, String} =
        config["impurity_density_profile"]["annealing_temperature_ranges"]["high_T"]["D0"],
    highT_H::Union{RealQuantity, String} =
        config["impurity_density_profile"]["annealing_temperature_ranges"]["high_T"]["H"],
    a::Union{Real, String} =
        config["experimental_parameters"]["a"],
    b::Union{Real, String} =
        config["experimental_parameters"]["b"],
    input_units::Union{Missing, NamedTuple} = missing
)

    temperature_unit = !ismissing(input_units) ? input_units.temperature : u"K"
    diffusivity_unit = !ismissing(input_units) ?
        input_units.length^2 / internal_time_unit : u"m^2/s"

    lowT = LithiumDiffusionBranch{T}(
        _parse_value(T, lowT_Tmin, temperature_unit),
        _parse_value(T, lowT_Tmax, temperature_unit),
        _parse_value(T, lowT_D0,  diffusivity_unit),
        _parse_value(T, lowT_H,   u"cal")
    )

    highT = LithiumDiffusionBranch{T}(
        _parse_value(T, highT_Tmin, temperature_unit),
        _parse_value(T, highT_Tmax, temperature_unit),
        _parse_value(T, highT_D0,  diffusivity_unit),
        _parse_value(T, highT_H,   u"cal")
    )

    diffusion_params = LithiumDiffusionParameters{T}(lowT, highT)

    saturation_params = LithiumSaturationParameters{T}(
        T(a),
        T(b)
    )

    return ThermalDiffusionLithiumDensityParameters{T}(
        diffusion_params,
        saturation_params
    )
end


function ThermalDiffusionLithiumParameters(config::AbstractDict, input_units::Union{Missing, NamedTuple} = missing; kwargs...)
    if !haskey(config, "impurity_density_profile") 
        throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'impurity_density_profile'.")) 
    elseif !haskey(config["impurity_density_profile"], "annealing_temperature_ranges") 
        throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'impurity_density_profile/annealing_temperature_ranges'.")) 
    elseif !haskey(config["impurity_density_profile"]["annealing_temperature_ranges"], "low_T")
        throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'impurity_density_profile/annealing_temperature_ranges/low_T'."))
    elseif !haskey(config["impurity_density_profile"]["annealing_temperature_ranges"], "high_T")
        throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'impurity_density_profile/annealing_temperature_ranges/high_T'."))
    end

    for axis in ("T_min", "T_max", "D0", "H")
        if !haskey(config["impurity_density_profile"]["annealing_temperature_ranges"]["low_T"], axis)
            throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'impurity_density_profile/annealing_temperature_ranges/low_T/$(axis)'."))
        end
        if !haskey(config["impurity_density_profile"]["annealing_temperature_ranges"]["high_T"], axis)
            throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'impurity_density_profile/annealing_temperature_ranges/high_T/$(axis)'."))
        end
    end

    if !haskey(config, "experimental_parameters") 
        throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'experimental_parameters'.")) 
    elseif !haskey(config["experimental_parameters"], "a") 
        throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'experimental_parameters/a'.")) 
    elseif !haskey(config["experimental_parameters"], "b") 
        throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'experimental_parameters/b'.")) 
    end

    _ThermalDiffusionLithiumParameters(config; input_units = input_units, kwargs...)
end

const default_ThermalDiffusionLithiumDensity_config_file = joinpath(get_path_to_example_config_files(), "ThermalDiffusionLithiumDensity/Thermal_Diffusion_Config.yaml")
function ThermalDiffusionLithiumParameters(config_filename::AbstractString = default_ThermalDiffusionLithiumDensity_config_file, input_units::Union{Missing, NamedTuple} = missing; kwargs...)
    ThermalDiffusionLithiumParameters(parse_config_file(config_filename), input_units; kwargs...)
end