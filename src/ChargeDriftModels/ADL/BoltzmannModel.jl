mutable struct BoltzmannModel{T <: SSDFloat} <: AbstractTemperatureModel{T}

    #temperatures needed for the drift velocity rescaling
    temperature::T      # actual temperature of the crystal
    reftemperature::T   # temperature at which the reference values were measured

    #fit parameters for electrons (e) and holes (h) for <100> and <111> axes
    p1e100::T
    p2e100::T
    p3e100::T

    p1e111::T
    p2e111::T
    p3e111::T

    p1h100::T
    p2h100::T
    p3h100::T

    p1h111::T
    p2h111::T
    p3h111::T
end

function BoltzmannModel(config_file::Dict; T::Type{<:AbstractFloat} = Float32)::BoltzmannModel
    return BoltzmannModel{T}(config_file::Dict)
end

function BoltzmannModel{T}(config_file::Dict; temperature::Union{Missing, T} = missing)::BoltzmannModel where T <: SSDFloat
    config_file["temperature_dependence"]["model"] != "Boltzmann" ? error() : nothing
    ismissing(temperature) ? temperature = config_file["temperature_dependence"]["temperature"] : nothing
    m = BoltzmannModel{T}(
        temperature, #temperature = config_file["temperature_dependence"]["temperature"]
        config_file["drift"]["velocity"]["temperature"],

        config_file["temperature_dependence"]["parameters"]["e100"]["p1"],
        config_file["temperature_dependence"]["parameters"]["e100"]["p2"],
        config_file["temperature_dependence"]["parameters"]["e100"]["p3"],

        config_file["temperature_dependence"]["parameters"]["e111"]["p1"],
        config_file["temperature_dependence"]["parameters"]["e111"]["p2"],
        config_file["temperature_dependence"]["parameters"]["e111"]["p3"],

        config_file["temperature_dependence"]["parameters"]["h100"]["p1"],
        config_file["temperature_dependence"]["parameters"]["h100"]["p2"],
        config_file["temperature_dependence"]["parameters"]["h100"]["p3"],

        config_file["temperature_dependence"]["parameters"]["h111"]["p1"],
        config_file["temperature_dependence"]["parameters"]["h111"]["p2"],
        config_file["temperature_dependence"]["parameters"]["h111"]["p3"]
    )
end

function scale_to_given_temperature(m::BoltzmannModel{T})::NTuple{4,T} where T <: SSDFloat

    function f(p1::T, p2::T, p3::T)::T where T <: SSDFloat
        return (p1+ p2*exp(-p3/m.reftemperature))/(p1 + p2*exp(-p3/m.temperature))
    end

    scale_e100 = f(m.p1e100, m.p2e100, m.p3e100)
    scale_e111 = f(m.p1e111, m.p2e111, m.p3e111)
    scale_h100 = f(m.p1h100, m.p2h100, m.p3h100)
    scale_h111 = f(m.p1h111, m.p2h111, m.p3h111)

    return scale_e100, scale_e111, scale_h100, scale_h111
end
