@doc raw"""

    struct PowerLawTemperatureModel{T <: SSDFloat}

Values needed to parametrize the temperature dependence of the longitudinal drift velocity 
of electrons or holes along a crystal axis as a function of the electric field strength.
This model assumes a power law dependence of the mobility with temperature.

## Background information 

The parametrization for the temperature dependence of drift velocity, as a function of the electric 
field strength, $E$, was proposed by [M.A. Omar and L. Reggiani](https://www.sciencedirect.com/science/article/pii/0038110187900633):
```math
\quad\mu_0(T) = A/T^P~, \quad V_s(T) = B\tanh^{1/2}(\theta/2T) 
```
Note that the saturation velocity, $V_s$, is related to $E_0$ via
```math
E_0(T) = V_s(T)/\mu_0(T)
```
The four parameters, $A$, $P$, $B$, $\theta$, are different for electrons and holes. 
This model can be used to scale the $\mu_0$ and $E_0$ parameters of the [`VelocityParameters`](@ref) struct directly (assuming these where measured at a reference temperature of 77K) as follows:
```math
\mu_0(T) = \mu_0(77K) \left(\frac{T}{77K}\right)^{-P}~, \quad E_0(T) = E_0(77K) \sqrt{\frac{\tanh(\theta/2T)}{\tanh(\theta/2 \cdot 77K)}} \frac{\mu_0(77K)}{\mu_0(T)}
```

!!! note
    The work was confined to the fields along the <100> axis. 
    Here the same values are assumed for the <111> axis. 
    
## Fields
* `p_e::T`: Exponent $P$ for electrons.
* `theta_e::T`: Characteristic temperature $\theta$ for electrons.
* `p_h::T`: Exponent $P$ for holes.
* `theta_h::T`: Characteristic temperature $\theta$ for holes.

"""

mutable struct PowerLawTemperatureModel{T <: SSDFloat} <: AbstractTemperatureModel{T}
    reference_temperature::T
    p_e::T
    theta_e::T
    p_h::T
    theta_h::T
end

function _PowerLawTemperatureModel(
    config::AbstractDict; 
    T::Type = Float32,
    reference_temperature::Union{RealQuantity,String} = config["temperature_dependence"]["reference_temperature"],
    p_e::Union{Real,String} = config["temperature_dependence"]["parameters"]["e"]["p"],
    theta_e::Union{RealQuantity,String} = config["temperature_dependence"]["parameters"]["e"]["theta"],
    p_h::Union{Real,String} = config["temperature_dependence"]["parameters"]["h"]["p"],
    theta_h::Union{RealQuantity,String} = config["temperature_dependence"]["parameters"]["h"]["theta"]
    )::PowerLawTemperatureModel{T}

    reftemp = _parse_value(T, reference_temperature, u"K")
    p_e = _parse_value(T, p_e, NoUnits)
    theta_e = _parse_value(T, theta_e, u"K")
    p_h = _parse_value(T, p_h, NoUnits)
    theta_h = _parse_value(T, theta_h, u"K")

    return PowerLawTemperatureModel{T}(reftemp, p_e, theta_e, p_h, theta_h)
end

# Check the syntax of the DriftModel temperature dependence config file before parsing
function PowerLawTemperatureModel(config::AbstractDict; kwargs...)
    if !haskey(config, "temperature_dependence") 
         throw(ConfigFileError("PowerLawTemperatureModel config file needs entry 'temperature_dependence'."))
    elseif !haskey(config["temperature_dependence"], "reference_temperature") 
        throw(ConfigFileError("PowerLawTemperatureModel config file needs entry 'temperature_dependence/reference_temperature'.")) 
    elseif !haskey(config["temperature_dependence"], "parameters")
        throw(ConfigFileError("PowerLawTemperatureModel config file needs entry 'temperature_dependence/parameters'."))
    end

    for axis in ("e", "h")
        if !haskey(config["temperature_dependence"]["parameters"], axis)
            throw(ConfigFileError("PowerLawTemperatureModel config file needs entry 'temperature_dependence/parameters/$(axis)'."))
        else
            for param in ("p", "theta")
                if !haskey(config["temperature_dependence"]["parameters"][axis], param)
                    throw(ConfigFileError("PowerLawTemperatureModel config file needs entry 'temperature_dependence/parameters/$(axis)/$(param)'."))
                end
            end
        end
    end
     _PowerLawTemperatureModel(config; kwargs...)
end

function PowerLawTemperatureModel(config_filename::AbstractString = default_ADL2016_config_file; kwargs...)
    PowerLawTemperatureModel(parse_config_file(config_filename); kwargs...)
end

PowerLawTemperatureModel{T}(args...; kwargs...) where {T <: SSDFloat} = PowerLawTemperatureModel(args...; T=T, kwargs...)

# scale velocity parameters given PowerLawTemperatureModel parameters
function scale_to_temperature_powerlaw(v::VelocityParameters{T}, Temp::T, reftemp::T, p::T, theta::T) where {T <: SSDFloat}
    μ0::T = v.mu0 * (Temp / reftemp) ^ -p
    β::T  = v.beta
    E0::T = v.E0 * sqrt(tanh(theta / (2*Temp)) / tanh(theta/ (2 * reftemp))) * v.mu0 / μ0
    μn::T = μ0 / v.mu0 * v.mun
    return VelocityParameters{T}(μ0, β, E0, μn)
end

function scale_to_temperature(cdm::CDM, Temp::Union{<:Real,Unitful.Temperature}) where {T,M,N,CDM<:ADL2016ChargeDriftModel{T,M,N,PowerLawTemperatureModel{T}}}
    Temp = _parse_value(T, Temp, u"K")
    reftemp = cdm.temperaturemodel.reference_temperature
    theta_e = cdm.temperaturemodel.theta_e
    theta_h = cdm.temperaturemodel.theta_h
    p_e = cdm.temperaturemodel.p_e
    p_h = cdm.temperaturemodel.p_h
    # scale electron parameters
    electrons = scale_to_temperature_powerlaw(cdm.electrons, Temp, reftemp, p_e, theta_e)
    # scale hole parameters (take <100> values and assume no anisotropy)
    h100 = scale_to_temperature_powerlaw(cdm.holes.axis100, Temp, reftemp, p_h, theta_h)
    h111 = scale_to_temperature_powerlaw(cdm.holes.axis111, Temp, reftemp, p_h, theta_h)
    holes = CarrierParameters(h100, h111)
    CDM(electrons, holes, cdm.crystal_orientation, cdm.γ, cdm.parameters, cdm.temperaturemodel)
end

function scale_to_temperature(cdm::CDM, Temp::Union{<:Real,Unitful.Temperature}) where {T,M,N,CDM<:ADLChargeDriftModel{T,M,N,PowerLawTemperatureModel{T}}}
    Temp = _parse_value(T, Temp, u"K")
    reftemp = cdm.temperaturemodel.reference_temperature
    theta_e = cdm.temperaturemodel.theta_e
    theta_h = cdm.temperaturemodel.theta_h
    p_e = cdm.temperaturemodel.p_e
    p_h = cdm.temperaturemodel.p_h
    # scale electron parameters (take <100> values and assume no anisotropy)
    e100 = scale_to_temperature_powerlaw(cdm.electrons.axis100, Temp, reftemp, p_e, theta_e)
    e111 = scale_to_temperature_powerlaw(cdm.electrons.axis111, Temp, reftemp, p_e, theta_e)
    electrons = CarrierParameters(e100, e111)
    # scale hole parameters (take <100> values and assume no anisotropy)
    h100 = scale_to_temperature_powerlaw(cdm.holes.axis100, Temp, reftemp, p_h, theta_h)
    h111 = scale_to_temperature_powerlaw(cdm.holes.axis111, Temp, reftemp, p_h, theta_h)
    holes = CarrierParameters(h100, h111)
    CDM(electrons, holes, cdm.crystal_orientation, cdm.γ, cdm.parameters, cdm.temperaturemodel)
end
