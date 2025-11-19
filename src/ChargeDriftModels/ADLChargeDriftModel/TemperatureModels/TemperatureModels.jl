include("PowerLawTemperatureModel.jl")
include("VacuumTemperatureModel.jl")

scale_to_temperature(cdm::AbstractChargeDriftModel, ::Missing) = cdm

@inline function TemperatureModel(T::DataType, config::AbstractDict, input_units::Union{Missing, NamedTuple} = missing)
    model::String = if haskey(config, "temperature_dependence") && haskey(config["temperature_dependence"], "model")
        _model::String = config["temperature_dependence"]["model"] 
        if _model in ("Omar1987", "Vacuum")
            _model
        elseif _model in ("Linear", "Boltzmann", "PowerLaw")
            @warn "Since v0.11.0, temperature scaling is only supported at the drift parameter level. The now deprecated models, \"Linear\", \"Boltzmann\" and \"PowerLaw\", relied on the scaling of longitudinal drift velocity at every step of charge carrier drift. This scaling is no longer supported. Please use \"Omar1987\" for the `PowerLawTemperatureModel`, which directly scales charge drift parameters. \nThe drift parameters will not be rescaled."
            "Vacuum"
        else
            @info "Config File does not suit any of the predefined temperature models. The drift parameters will not be rescaled."
            "Vacuum"
        end
    else
        "Vacuum"
    end
    TemperatureModel(T, Val{Symbol(model)}(), config, input_units)
end