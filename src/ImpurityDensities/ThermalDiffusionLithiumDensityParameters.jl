"""
    struct LithiumDiffusionParameters{T <: SSDFloat}

Set of parameters to calculate the lithium diffusivity `D` at a given annealing temperature `T_an` using
```math
D(T_\text{an}) = D_0 \exp(-H / RT+\text{an})
```

## Parametric types
* `T`: Precision type.
    
## Fields
* `T_min::T`: Minimum temperature (in K) for which the parameters are valid.
* `T_max::T`: Maximum temperature (in K) for which the parameters are valid.
* `D0::T`: Diffusivity constant in m^2/s.
* `H::T`: Activation energy in cal/mol.
"""
struct LithiumDiffusionParameters{T <: SSDFloat}
    T_min::T
    T_max::T
    D0::T
    H::T
end

function calculate_lithium_diffusivity(lithium_annealing_temperature::T, parameters::NTuple{N, LithiumDiffusionParameters{T}})::T where {T <: SSDFloat, N}
    if !(parameters[1].T_min ≤ lithium_annealing_temperature ≤ parameters[end].T_max)
        throw(ArgumentError("Invalid lithium_annealing_temperature=$(lithium_annealing_temperature): expected $(parameters[1].T_min) ≤ lithium_annealing_temperature ≤ $(parameters[end].T_max)."))
    end

    for parameter in parameters
        if lithium_annealing_temperature <= parameter.T_max
            return parameter.D0 * exp(-parameter.H/(R_gas*lithium_annealing_temperature))
        end
    end
end


"""
    struct LithiumSaturationParameters{T <: SSDFloat}

Set of experimental fit parameters to calculate the lithium saturated density `N_s` in 1/cm³ at a given annealing temperature `T_an` using
```math
N_s(T_\text{an}) = 10^{a - b/T_\text{an})
```

## Parametric types
* `T`: Precision type.
    
## Fields
* `a::T`: Fit parameter.
* `b::T`: Fit parameter.`
"""
struct LithiumSaturationParameters{T <: SSDFloat}
    a::T
    b::T
end

function calculate_lithium_saturated_density(lithium_annealing_temperature::T, parameters::LithiumSaturationParameters{T})::T where {T <: SSDFloat}
    ustrip(internal_length_unit^-3, exp10(parameters.a - parameters.b/lithium_annealing_temperature) * u"cm^-3")
end


"""
    ThermalDiffusionLithiumDensityParameters{T <: SSDFloat,N}

Set of parameters to calculate the lithium diffusivity and the lithium saturated density.

## Parametric types
* `T`: Precision type.
* `N`: Number of temperature ranges to calculate the lithium diffusivity.
    
## Fields
* `diffusion::NTuple{N, LithiumDiffusionParameters{T}}`: Lithium diffusivity parameters.
* `saturation::LithiumSaturationParameters{T}`: Lithium saturated density parameters.`

See also [`LithiumDiffusionParameters`](@ref) and [`LithiumSaturationParameters`](@ref).
"""
struct ThermalDiffusionLithiumDensityParameters{T <: SSDFloat,N}
    diffusion::NTuple{N, LithiumDiffusionParameters{T}}
    saturation::LithiumSaturationParameters{T}
end


function _ThermalDiffusionLithiumParameters(
        config::AbstractDict;
        T::Type = Float32,
        temperature_ranges::AbstractVector{<:AbstractDict} = config["annealing_temperature_ranges"] ,
        a::Union{Real, String} = config["experimental_parameters"]["a"],
        b::Union{Real, String} = config["experimental_parameters"]["b"],
        input_units::Union{Missing, NamedTuple} = missing
    )

    temperature_unit = !ismissing(input_units) ? input_units.temperature : internal_temperature_unit
    diffusivity_unit = !ismissing(input_units) ? input_units.length^2 / internal_time_unit : internal_length_unit^2/internal_time_unit

    diffusion_params = Tuple(
        LithiumDiffusionParameters{T}(
            _parse_value(T, r["T_min"], temperature_unit),
            _parse_value(T, r["T_max"], temperature_unit),
            _parse_value(T, r["D0"], diffusivity_unit),
            _parse_value(T, r["H"], internal_activation_energy_unit)
        )
        for r in temperature_ranges
    )

    for i in eachindex(diffusion_params)
        if diffusion_params[i].T_max < diffusion_params[i].T_min
            throw(ConfigFileError("Invalid annealing temperature range $i: T_max must be ≥ T_min."))
        end
        if i > 1 && diffusion_params[i-1].T_max != diffusion_params[i].T_min
            throw(ConfigFileError("Annealing temperature ranges must be contiguous and increasing: range $(i-1) ends at a different temperature than range $i starts."))
        end
    end

    saturation_params = LithiumSaturationParameters{T}(_parse_value(T, a, Unitful.NoUnits), _parse_value(T, b, Unitful.NoUnits))
    
    ThermalDiffusionLithiumDensityParameters(diffusion_params, saturation_params)
end


function ThermalDiffusionLithiumParameters(config::AbstractDict, input_units::Union{Missing, NamedTuple} = missing; kwargs...)
    if !haskey(config, "annealing_temperature_ranges") 
        throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'annealing_temperature_ranges'.")) 
    elseif isempty(config["annealing_temperature_ranges"])
        throw(ConfigFileError("ThermalDiffusionLithiumDensity config file cannot have an empty 'annealing_temperature_ranges'."))
    end 

    for (i, temperature_range) in enumerate(config["annealing_temperature_ranges"])
        for axis in ("T_min", "T_max", "D0", "H")
            if !haskey(temperature_range, axis)
                throw(ConfigFileError("ThermalDiffusionLithiumDensity config file needs entry 'annealing_temperature_ranges[$i]/$(axis)'."))
            end
        end
        for k in keys(temperature_range)
            if !(k in ("T_min", "T_max", "D0", "H"))
                @warn "Encountered unexpected key `$(k)` in 'annealing_temperature_ranges[$i]': Allowed keys are `T_min`, `T_max`, `D0` and `H`. Key `$(k)` will be ignored."
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