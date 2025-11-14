@doc raw"""

    struct PowerLawModel{T <: SSDFloat}

Values needed to parametrize the temperature dependence of the longitudinal drift velocity 
of electrons or holes along a crystal axis as a function of the electric field strength.
This model assumes a power law dependence of the mobility with temperature.

## Background information 

The parametrization for the temperature dependence of drift velocity, as a function of the electric 
field strength, ``E``, was proposed by [M.A. Omar and L. Reggiani](https://www.sciencedirect.com/science/article/pii/0038110187900633):
```math
\quad\mu_0(T) = A/T^P~, \quad V_s(T) = B\tanh^{1/2}(\theta/2T) 
```
Note that the saturation velocity, ``V_s``, is related to ``E_0`` via
```math
E_0(T) = V_s(T)/\mu_0(T)
```
The four parameters, ``A``, ``P``, ``B``, ``\theta``, are different for electrons and holes. 
This model can be used to scale the ``\mu_0`` and ``E_0`` parameters of the [`VelocityParameters`](@ref) struct directly (assuming these where measured at a reference temperature of 77K) as follows:
```math
\mu_0(T) = \mu_0(77K) \left(\frac{T}{77K}\right)^{-P}~, \quad E_0(T) = E_0(77K) \sqrt{\frac{\tanh(\theta/2T)}{\tanh(\theta/2 \cdot 77K)}} \frac{\mu_0(77K)}{\mu_0(T)}
```

!!! note
    The work was confined to the fields along the <100> axis. 
    Here the same values are assumed for the <111> axis. 
    
## Fields
* `p_e100::T`: Exponent ``P`` for electrons along the <100> axis.
* `theta_e100::T`: Characteristic temperature ``\theta`` for electrons along the <100> axis.
* `p_h100::T`: Exponent ``P`` for holes along the <100> axis.
* `theta_h100::T`: Characteristic temperature ``\theta`` for holes along the <100> axis.

"""

mutable struct PowerLawModel{T <: SSDFloat} <: AbstractTemperatureModel{T}
    reference_temperature::T
    p_e100::T
    theta_e100::T
    p_h100::T
    theta_h100::T
end

function PowerLawModel(config_file::AbstractDict; T::Type{<:AbstractFloat} = Float32)::PowerLawModel
    return PowerLawModel{T}(config_file::AbstractDict)
end

function PowerLawModel{T}(config_file::AbstractDict)::PowerLawModel where {T <: SSDFloat}
    config_file["temperature_dependence"]["model"] != "PowerLaw" ? error() : nothing
    m = PowerLawModel{T}(
        _parse_value(T, config_file["temperature_dependence"]["reference_temperature"], u"K"),
        _parse_value(T, config_file["temperature_dependence"]["parameters"]["e100"]["p"], NoUnits),
        _parse_value(T, config_file["temperature_dependence"]["parameters"]["e100"]["theta"], u"K"),
        _parse_value(T, config_file["temperature_dependence"]["parameters"]["h100"]["p"], NoUnits),
        _parse_value(T, config_file["temperature_dependence"]["parameters"]["h100"]["theta"], u"K")
    )
end

function scale_to_temperature(cdm::ADL2016ChargeDriftModel{T,M,N,PowerLawModel{T}}, Temp::RealQuantity) where {T,M,N}
    Temp = _parse_value(T, Temp, u"K")
    reftemp = cdm.temperaturemodel.reference_temperature
    theta_e = cdm.temperaturemodel.theta_e100
    theta_h = cdm.temperaturemodel.theta_h100
    p_e = cdm.temperaturemodel.p_e100
    p_h = cdm.temperaturemodel.p_h100

    # scale electron parameters
    μe::T = cdm.electrons.mu0 * (Temp / reftemp) ^ -p_e
    βe::T = cdm.electrons.beta
    Ee::T = cdm.electrons.E0 * sqrt(tanh(theta_e / (2*Temp)) / tanh(theta_e/ (2 * reftemp))) * cdm.electrons.mu0 / μe
    μn::T = μe / cdm.electrons.mu0 * cdm.electrons.mun
    electrons = SolidStateDetectors.VelocityParameters{T}(μe, βe, Ee, μn)
    
    # scale hole parameters (take <100> values and assume no anisotropy for now)
    μh::T = cdm.holes.axis100.mu0 * (Temp / reftemp) ^ -p_h
    βh::T = cdm.holes.axis100.beta
    Eh::T = cdm.holes.axis100.E0 * sqrt(tanh(theta_h / (2*Temp)) / tanh(theta_h / (2 * reftemp))) * cdm.holes.axis100.mu0 / μh
    h100 = SolidStateDetectors.VelocityParameters{T}(μh, βh, Eh, zero(T))
    μh = cdm.holes.axis111.mu0 * (Temp / reftemp) ^ -p_h
    βh = cdm.holes.axis111.beta
    Eh = cdm.holes.axis111.E0 * sqrt(tanh(theta_h / (2*Temp)) / tanh(theta_h / (2 * reftemp))) * cdm.holes.axis111.mu0 / μh
    h111 = SolidStateDetectors.VelocityParameters{T}(μh, βh, Eh, zero(T))
    holes = SolidStateDetectors.CarrierParameters(h100, h111)
    
    ADL2016ChargeDriftModel{T,M,N,PowerLawModel{T}}(electrons, holes, cdm.crystal_orientation, cdm.γ, cdm.parameters, cdm.temperaturemodel) 
end  
