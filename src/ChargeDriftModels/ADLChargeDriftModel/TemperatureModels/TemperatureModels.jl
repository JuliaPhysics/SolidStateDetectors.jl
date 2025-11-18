include("PowerLawTemperatureModel.jl")
include("VacuumTemperatureModel.jl")

scale_to_temperature(cdm::AbstractChargeDriftModel, ::Missing) = cdm

@inline function TemperatureModel(T::DataType, config::AbstractDict)
     model::String = if haskey(config, "temperature_dependence") && haskey(config["temperature_dependence"], "model")
        config["temperature_dependence"]["model"]
    else
        "Vacuum"
    end
    model_set::String = if model in ("Omar1987", "Vacuum") 
        model
    else 
        if model in ("Linear", "Boltzmann", "PowerLaw")
            @warn "Since v0.11.0, temperature scaling is only supported at the drift parameter level. The now deprecated models, \"Linear\", \"Boltzmann\" and \"PowerLaw\", relied on the scaling of longitudinal drift velocity at every step of charge carrier drift. This scaling is no longer supported. Please use \"Omar1987\" for the `PowerLawTemperatureModel`, which directly scales charge drift parameters. \nThe drift parameters will not be rescaled."
        else
            @info "Config File does not suit any of the predefined temperature models. The drift parameters will not be rescaled."
        end
        "Vacuum"
    end
    TemperatureModel(T, Val{Symbol(model_set)}(), config)
end