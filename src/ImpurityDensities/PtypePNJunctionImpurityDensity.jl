"""
struct ThermalDiffusionLithiumDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

A PN junction impurity model based on lithium thermal diffusion and constant bulk impurity density.
The surface lithium density is saturated.

ref: [Dai _et al._ (2023)](https://doi.org/10.1016/j.apradiso.2022.110638)
 
## Fields
* `surface_imp_model::ThermalDiffusionLithiumDensity{T,G}`: the density profile of lithium (n-type).
* `bulk_imp_model::AbstractImpurityDensity{T}`: the density of the p-type impurity.
"""

struct PtypePNJunctionImpurityDensity{T <: SSDFloat, G <: Union{<:AbstractGeometry, Nothing}} <: AbstractImpurityDensity{T}
    surface_imp_model::ThermalDiffusionLithiumDensity{T, G}
    bulk_imp_model::AbstractImpurityDensity{T}
end
function PtypePNJunctionImpurityDensity{T}(
    lithium_annealing_temperature::T,
    lithium_annealing_time::T,
    contact_with_lithium_doped::G,
    bulk_imp_model::AbstractImpurityDensity{T},
    distance_to_contact::Function = pt::AbstractCoordinatePoint{T} -> ConstructiveSolidGeometry.distance_to_surface(pt, contact_with_lithium_doped)
    ) where {T <: SSDFloat, G <: Union{<:AbstractGeometry, Nothing}}
    surface_imp_model = ThermalDiffusionLithiumDensity{T}(lithium_annealing_temperature, lithium_annealing_time, contact_with_lithium_doped, distance_to_contact=distance_to_contact)
    PtypePNJunctionImpurityDensity{T,G}(surface_imp_model, bulk_imp_model)
end

function ImpurityDensity(T::DataType, t::Val{:PtypePNjunction}, dict::AbstractDict, input_units::NamedTuple)
    lithium_annealing_temperature = _parse_value(T, get(dict, "lithium_annealing_temperature", 623u"K"), input_units.temperature)
    lithium_annealing_time = _parse_value(T, get(dict, "lithium_annealing_time", 18u"minute"), internal_time_unit)
    contact_with_lithium_doped = get(dict, "contact_with_lithium_doped", nothing)
    bulk_imp_model = haskey(dict, "bulk_imp_model") ? ImpurityDensity(T, dict["bulk_imp_model"], input_units) : ConstantImpurityDensity{T}(-1e16)
    PtypePNJunctionImpurityDensity{T}(lithium_annealing_temperature, lithium_annealing_time, contact_with_lithium_doped, bulk_imp_model)
end

function get_impurity_density(PtypePNjunction::PtypePNJunctionImpurityDensity{T, G}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat, G <: Union{<:AbstractGeometry, Nothing}}
    get_impurity_density(PtypePNjunction.bulk_imp_model, pt)+get_impurity_density(PtypePNjunction.surface_imp_model, pt)
end
