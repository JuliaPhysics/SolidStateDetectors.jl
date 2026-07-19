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

# Hole mobility: Dai et al. (2023), eqs. 5-8, converted from their cm-based
# units to SI (densities m^-3, mobilities m^2/V/s). The neutral-impurity term
# keeps the Sclar form 0.82 * (light+heavy hole Erginsoy constants) * f(T); it
# collapses to 4.455e21/Nn * (T^0.5 + 4.281*T^-0.5), i.e. 4.46e19 in the
# paper's cm-units - the paper's printed 4.46e29 is an exponent typo (with it,
# no physical Nn could reproduce the mobility matching described there).
function _calculate_mobility_with_impurities(
    Nn::T, bulk_imp::T, surface_imp::T,
    temperature::T,
    ::Type{Hole})::T where {T}
    Ni::T = abs(-bulk_imp + surface_imp)

    μI::T = 2.35e19*temperature^1.5/Ni/log(9.13e19*temperature^2/Ni) + 1.51e20*temperature^1.5/Ni/log(5.82e20*temperature^2/Ni)
    μA::T = 7.77e3 * temperature^-1.5
    μN::T = 1e2/Nn * (2.31e18+2.36e20) * 0.82 * (0.228*temperature^0.5 + 0.976*temperature^-0.5)

    1/(1/μI + 1/μA + 1/μN)
end

# Electron mobility: Ma et al. (arXiv:1705.09562), Erginsoy (1.07e20 cm-units
# -> 1.07e22 SI) with Sclar's temperature correction (0.28*T^0.5 + 0.54*T^-0.5)
# and Bardeen-Shockley acoustic term (9.32e7 cm-units -> 9.32e3 SI).
function _calculate_mobility_with_impurities(
    Nn::T, bulk_imp::T, surface_imp::T,
    temperature::T,
    ::Type{Electron})::T where {T}
    Ni::T = abs(-bulk_imp + surface_imp)

    μI::T = 2.442e20*temperature^1.5/Ni/(log(2.496e20*temperature^2/Ni))
    μA::T = 9.32e3 * temperature^-1.5
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
