mutable struct PowerLawModel{T <: SSDFloat} <: AbstractTemperatureModel{T}

    #temperatures needed for the drift velocity rescaling
    temperature::T      # actual temperature of the crystal
    reftemperature::T   # temperature at which the reference values were measured

    #fit parameters for electrons (e) and holes (h) for <100> and <111> axes
    p1e100::T
    p1e111::T
    p1h100::T
    p1h111::T
end

function PowerLawModel(config_file::Dict; T::Type{<:AbstractFloat} = Float32)::PowerLawModel
    return PowerLawModel{T}(config_file::Dict)
end

function PowerLawModel{T}(config_file::Dict; temperature::Union{Missing, T} = missing)::PowerLawModel where {T <: SSDFloat}
    config_file["temperature_dependence"]["model"] != "PowerLaw" ? error() : nothing
    ismissing(temperature) ? temperature = config_file["temperature_dependence"]["temperature"] : nothing
    m = PowerLawModel{T}(
        temperature, #temperature = config_file["temperature_dependence"]["temperature"]
        config_file["drift"]["velocity"]["temperature"],
        config_file["temperature_dependence"]["parameters"]["e100"]["p1"],
        config_file["temperature_dependence"]["parameters"]["e111"]["p1"],
        config_file["temperature_dependence"]["parameters"]["h100"]["p1"],
        config_file["temperature_dependence"]["parameters"]["h111"]["p1"]
    )
end

function scale_to_given_temperature(m::PowerLawModel{T})::NTuple{4,T} where {T <: SSDFloat}

    function f(p1::T = 0.0)::T where {T <: SSDFloat}
        return (m.reftemperature/m.temperature)^1.5
    end

    scale_e100 = f(m.p1e100)
    scale_e111 = f(m.p1e111)
    scale_h100 = f(m.p1h100)
    scale_h111 = f(m.p1h111)

    return scale_e100, scale_e111, scale_h100, scale_h111
end
