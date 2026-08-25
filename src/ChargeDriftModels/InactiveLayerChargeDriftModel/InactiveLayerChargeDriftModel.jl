"""
    struct InactiveLayerChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
        
Charge drift model in which the electrons and holes drift along the electric field.
Three factors are considered in the mobility calculation: ionized impurities, neutral impurities, and acoustic phonons.

Ref: [Dai _et al._ (2023)](https://doi.org/10.1016/j.apradiso.2022.110638)

## Fields
- `temperature::T`: the crystal temperature (Kelvin), the default is 90 (Kelvin).
- `neutral_imp_model::AbstractImpurityDensity{T}`: the neutral impurity density model. The default is a constant impurity density of 1e21 m⁻³.
- `bulk_imp_model::AbstractImpurityDensity{T}`: the bulk impurity density model. The default is the defined (bulk) impurity density if the model is constructed in the config, otherwise the default is a constant impurity density of -1e16 m⁻³.
- `surface_imp_model::AbstractImpurityDensity{T}`: the surface impurity density model. The default is the defined surface impurity density if the model is constructed in the config, otherwise the default is a constant impurity density of 0 m⁻³.
"""
struct InactiveLayerChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
    temperature::T
    neutral_imp_model::AbstractImpurityDensity{T}
    bulk_imp_model::AbstractImpurityDensity{T}
    surface_imp_model::AbstractImpurityDensity{T}
    
    # constructor with default values
    InactiveLayerChargeDriftModel{T}(
        temperature::T = T(90),
        neutral_imp_model::AbstractImpurityDensity{T} = ConstantImpurityDensity{T}(T(1e21)),
        bulk_imp_model::AbstractImpurityDensity{T} = ConstantImpurityDensity{T}(T(-1e16)),
        surface_imp_model::AbstractImpurityDensity{T} = ConstantImpurityDensity{T}(T(0))
    ) where {T <: SSDFloat} = new{T}(temperature, neutral_imp_model, bulk_imp_model, surface_imp_model)
end


function calculate_mobility(cdm::InactiveLayerChargeDriftModel{T}, pt::AbstractCoordinatePoint{T}, CC::Type{CC_type}) where {T <: SSDFloat, CC_type <: ChargeCarrier}
    _calculate_mobility_with_impurities(
        get_impurity_density(cdm.neutral_imp_model, pt),
        get_impurity_density(cdm.bulk_imp_model, pt),
        get_impurity_density(cdm.surface_imp_model, pt),
        cdm.temperature, CC)
end

function _calculate_mobility_with_impurities(
    Nn::T, bulk_imp::T, surface_imp::T,
    temperature::T,
    ::Type{Hole})::T where {T}
    Ni::T = abs(-bulk_imp + surface_imp)

    # Based on Dai et al. (2023) https://doi.org/10.1016/j.apradiso.2022.110638
    # - Eq. (6): But adjusting the exponents to account for conversion from cm <-> m
    μI::T = 2.35e19*temperature^1.5/Ni/log(9.13e19*temperature^2/Ni) + 1.51e20*temperature^1.5/Ni/log(5.82e20*temperature^2/Ni)

    # Based on Mei et al. (2016) https://doi.org/10.1088/1748-0221/11/12/P12021
    # - µA = 1/17.05 * 1.15x10^9 + 16.05/17.05 * 1.12x10^7 = 7.77x10^7
    # - Dividing the number by 10^4 to convert from cm^2/Vs to m^2/Vs (SI units)
    μA::T = 7.77e3 * temperature^-1.5

    # Based on Mei et al. (2016) https://doi.org/10.1088/1748-0221/11/12/P12021 and references therein:
    # - Eqs. (2.22)-(2.24): Combining light and heavy holes: 
    #   1/17.05 * 3.94x10^{19} + 16.05/17.05 * 2.51x10^{20} = 2.31x10^18 + 2.36x10^20
    # - Eq. (2.25): Multiplying with 0.82, taken from Eq.(6) in McGill & Baron (1974) https://doi.org/10.1103/PhysRevB.11.5208 
    # - Eq. (2.25): Temperature dependence as shown in Eq.(19) in Sclar (1970) https://doi.org/10.1103/PhysRev.104.1559
    #   assuming E_N = 0.71eV * β/D^2 with β = 1/17.05 * 0.044 + 16.05/17.05 * 0.28 = 0.266 (ratio between effective mass and free mass), 
    #   and D = 16.0 (relative permittivity) => E_N = 0.71eV * 0.266 / 16.0^2 = 0.7377meV
    #   The two terms are ustrip(u"K^(-1/2)", 2/3 * sqrt(Unitful.k / E_N)) = 0.22778018408674738
    #   and ustrip(u"K^(1/2)", 1/3 * sqrt(EN / Unitful.k)) = 0.9755994495886063
    # - The first factor 1e2 makes sure to convert the equation to SI units 
    #   (1e6 for converting Nn from cm^-3 to m^-3, and 1e-4 to convert µN from cm^2/Vs to m^2/Vs)
    μN::T = 1e2/Nn * (2.31e18+2.36e20) * 0.82 * (0.228*temperature^0.5 + 0.976*temperature^-0.5)

    1/(1/μI + 1/μA + 1/μN)
end

function _calculate_mobility_with_impurities(
    Nn::T, bulk_imp::T, surface_imp::T,
    temperature::T,
    ::Type{Electron})::T where {T}
    Ni::T = abs(-bulk_imp + surface_imp)

    # Based on Zhang et al. (2026) https://doi.org/10.1140/epjc/s10052-026-15508-3
    # - Eq. (6): But adjusting the exponents to account for conversion from cm <-> m
    μI::T = 2.442e20*temperature^1.5/Ni/(log(2.496e20*temperature^2/Ni))

    # Based on Mei et al. (2017) https://doi.org/10.1088/1748-0221/12/07/P07003
    # - Eq. (2.6): But adjusting the exponents to account for conversion from cm <-> m
    μA::T = 9.32e3 * temperature^-1.5

    # Based on Mei et al. (2017) https://doi.org/10.1088/1748-0221/12/07/P07003
    # - Eq. (2.7): First factor the temperature-independent value quoted as µN = 1.07x10^20 / Nn
    # - Eq. (2.8): Temperature dependence as shown in Eq.(19) in Sclar (1970) https://doi.org/10.1103/PhysRev.104.1559
    #   assuming E_N = 0.71eV *  β/D^2 with β = 0.12 (ratio between effective mass and free mass), and D = 16.0
    #   => E_N = 0.71eV * 0.12 / 16.0^2 = 0.333meV
    #   The two terms are 0.82 * ustrip(u"K^(-1/2)", 2/3 * sqrt(Unitful.k / E_N)) = 0.82 * 0.3392308735301575 = 0.27816931629472913
    #   and 0.82 * ustrip(u"K^(1/2)", 1/3 * sqrt(EN / Unitful.k)) = 0.82 * 0.6550766441441445 = 0.5371628481981984
    #   !! NOTE that the factor of 0.82 from Eq. (2.8) is already factored in with the coefficients 0.28 and 0.54
    μN::T = 1.07e22/Nn * (0.28*temperature^0.5 + 0.54*temperature^-0.5)

    1/(1/μI + 1/μA + 1/μN)
end

function InactiveLayerChargeDriftModel{T}(config::AbstractDict,
        imp_model::AbstractImpurityDensity, input_units::NamedTuple, temperature::RealQuantity
    ) where {T <: SSDFloat}

    temp::T = _parse_value(T, temperature, internal_temperature_unit)

    if haskey(config, "temperature")
        cdm_temperature::T = _parse_value(T, config["temperature"], input_units.temperature)
        if cdm_temperature != temp
            throw(ConfigFileError(
                "Temperature mismatch: InactiveLayerChargeDriftModel defines temperature = $(cdm_temperature) K, " *
                "but the semiconductor temperature is $(temp) K. " *
                "Define `temperature` in the semiconductor and either remove it from the charge drift model or make sure that the temperatures match."
            ))
        end
    end

    if temp < 50 || temp > 150
        @warn "Temperature = $(temp) K is outside the typical validated range (50–150 K)."
    end
    neutral_imp_model = if haskey(config, "neutral_impurity_density")
        ImpurityDensity(T, config["neutral_impurity_density"], input_units)
    else
        ConstantImpurityDensity{T}(1e21)
    end

    bulk_imp_model = if haskey(config, "bulk_impurity_density")
        ImpurityDensity(T, config["bulk_impurity_density"], input_units)
    elseif isdefined(imp_model, :bulk_imp_model)
        imp_model.bulk_imp_model
    else
        imp_model
    end

    surface_imp_model = if haskey(config, "surface_impurity_density")
        ImpurityDensity(T, config["surface_impurity_density"], input_units)
    elseif isdefined(imp_model, :surface_imp_model)
        imp_model.surface_imp_model
    else
        throw(ConfigFileError("There is no surface impurity density profile provided."))
    end
    InactiveLayerChargeDriftModel{T}(temp, neutral_imp_model, bulk_imp_model, surface_imp_model)
end

@fastmath function getVe(fv::SVector{3, T}, cdm::InactiveLayerChargeDriftModel{T}, current_pos::CartesianPoint{T} = zero(CartesianPoint{T})) where {T <: SSDFloat}
    -calculate_mobility(cdm, current_pos, Electron) * fv
end

@fastmath function getVh(fv::SVector{3, T}, cdm::InactiveLayerChargeDriftModel{T}, current_pos::CartesianPoint{T} = zero(CartesianPoint{T})) where {T <: SSDFloat}
    calculate_mobility(cdm, current_pos, Hole) * fv
end
