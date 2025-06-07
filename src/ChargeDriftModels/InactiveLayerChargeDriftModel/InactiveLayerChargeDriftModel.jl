"""
    struct InactiveLayerChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
        
Charge drift model in which the electrons and holes drift along the electric field.
There factors are considered in the mobility calculation: ionized impurities, neutral impurities, and acoustic phonon.

Ref: [Dai _et al._ (2023)](https://doi.org/10.1016/j.apradiso.2022.110638)

## Fields
- `calculate_mobility::Function`: Mobility calculation function
- `neutral_imp_model::AbstractImpurityDensity{T}`: the neutral impurity density model
- `bulk_imp_model::AbstractImpurityDensity{T}`: the bulk impurity density model (p-type impurity)
- `surface_imp_model::AbstractImpurityDensity{T}`: the surface impurity density model (n-type impurity)

## Extra field for constructing the model
- `temperature::T`: temperature of the crystal (Kelvin)
"""

struct InactiveLayerChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
    calculate_mobility::Function
	neutral_imp_model::AbstractImpurityDensity{T}
	bulk_imp_model::AbstractImpurityDensity{T}
    surface_imp_model::AbstractImpurityDensity{T}
end

function InactiveLayerChargeDriftModel{T}(
	temperature::T, # Kelvin
	neutral_imp_model::AbstractImpurityDensity{T},
	bulk_imp_model::AbstractImpurityDensity{T},
    surface_imp_model::AbstractImpurityDensity{T}
	) where {T <: SSDFloat}

	function calculate_mobility(pt::AbstractCoordinatePoint{T}, CC::Type{CC_type}) where {T <: SSDFloat, CC_type <: ChargeCarrier}
		calculate_mobility_(
			get_impurity_density(neutral_imp_model, pt),
			get_impurity_density(bulk_imp_model, pt),
			get_impurity_density(surface_imp_model, pt),
			temperature, CC)
	end

    InactiveLayerChargeDriftModel{T}(calculate_mobility, neutral_imp_model, bulk_imp_model, surface_imp_model)
end

function calculate_mobility_(
	Nn::T, bulk_imp::T, surface_imp::T,
	temperature::T,
	::Type{Hole}) where {T}
	Ni = -bulk_imp + surface_imp

	μI = 2.35e19*temperature^1.5/Ni/log(9.13e19*temperature^2/Ni) + 1.51e20*temperature^1.5/Ni/log(5.82e20*temperature^2/Ni)
	μA = 7.77e3 * temperature^-1.5
	μN = 1e2/Nn * (2.31e18+2.36e20) * 0.82 * (0.228*temperature^0.5 + 0.976*temperature^-0.5)

	1/(1/μI + 1/μA + 1/μN)
end
function calculate_mobility_(
	Nn::T, bulk_imp::T, surface_imp::T,
	temperature::T,
	::Type{Electron}) where {T}
	Ni = -bulk_imp + surface_imp

	μI = 2.442e20*temperature^1.5/Ni/(log(2.496e20*temperature^2/Ni))
	μA = 9.32e3 * temperature^-1.5
	μN = 1.07e22/Nn * (0.28*temperature^0.5 + 0.54*temperature^-0.5)

	1/(1/μI + 1/μA + 1/μN)
end

function InactiveLayerChargeDriftModel{T}(config::AbstractDict,
	imp_model::AbstractImpurityDensity, input_units::NamedTuple
	) where {T <: SSDFloat}
    temperature = _parse_value(T, get(config, "temperature", 90u"K"), input_units.temperature)

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
    InactiveLayerChargeDriftModel{T}(temperature, neutral_imp_model, bulk_imp_model, surface_imp_model)
end

@fastmath function getVe(fv::SVector{3, T}, cdm::InactiveLayerChargeDriftModel{T}, current_pos::CartesianPoint{T} = zero(CartesianPoint{T})) where {T <: SSDFloat}
    @inbounds begin 
        -cdm.calculate_mobility(current_pos, Electron)*fv
    end
end

@fastmath function getVh(fv::SVector{3, T}, cdm::InactiveLayerChargeDriftModel{T}, current_pos::CartesianPoint{T} = zero(CartesianPoint{T})) where {T <: SSDFloat}
    @inbounds begin 
        cdm.calculate_mobility(current_pos, Hole)*fv
    end
end