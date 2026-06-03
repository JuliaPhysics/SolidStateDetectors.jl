"""
struct PtypePNJunctionImpurityDensity{T <: SSDFloat, N} <: AbstractImpurityDensity{T}

A PN junction impurity model based on lithium thermal diffusion and constant bulk impurity density.
The surface lithium density is at the saturation level. Ref: [Dai _et al._ (2023)](https://doi.org/10.1016/j.apradiso.2022.110638)

## Parametric types 
* `T`: Precision type.
* `N`: Number of lithium-drifted contacts.
 
## Fields
* `surface_imp_model::ThermalDiffusionLithiumDensity{T,N}`: the density profile of lithium (n-type). Detailed in [`ThermalDiffusionLithiumDensity`](@ref).
* `bulk_imp_model::AbstractImpurityDensity{T}`: the density profile of the p-type impurities.
"""
struct PtypePNJunctionImpurityDensity{T <: SSDFloat, N} <: AbstractImpurityDensity{T}
    surface_imp_model::ThermalDiffusionLithiumDensity{T,N}
    bulk_imp_model::AbstractImpurityDensity{T}
    PtypePNJunctionImpurityDensity{T}(
        surface_imp_model::ThermalDiffusionLithiumDensity{T,N}, 
        bulk_imp_model::AbstractImpurityDensity{T}
    ) where {T,N} = new{T,N}(surface_imp_model, bulk_imp_model)
end
function PtypePNJunctionImpurityDensity{T}(
    lithium_annealing_temperature::T,
    lithium_annealing_time::T,
    contact_with_lithium_doped::G,
    inactive_contact_id::Union{<:Integer, <:NTuple{<:Any,<:Integer}, <:AbstractVector{<:Integer}},
    bulk_imp_model::AbstractImpurityDensity{T},
    distance_to_contact::Function = pt::AbstractCoordinatePoint{T} -> ConstructiveSolidGeometry.distance_to_surface(pt, contact_with_lithium_doped)
    ) where {T <: SSDFloat, G <: Union{<:AbstractGeometry, Nothing}}
    surface_imp_model = ThermalDiffusionLithiumDensity{T}(lithium_annealing_temperature, lithium_annealing_time, contact_with_lithium_doped, inactive_contact_id, distance_to_contact=distance_to_contact)
    PtypePNJunctionImpurityDensity{T}(surface_imp_model, bulk_imp_model)
end

function ImpurityDensity(T::DataType, t::Val{:PtypePNjunction}, dict::AbstractDict, input_units::NamedTuple)
    lithium_annealing_temperature = _parse_value(T, get(dict, "lithium_annealing_temperature", 623u"K"), input_units.temperature)
    lithium_annealing_time = _parse_value(T, get(dict, "lithium_annealing_time", 18u"minute"), internal_time_unit)
    contact_with_lithium_doped = get(dict, "contact_with_lithium_doped", nothing)
    inactive_contact_id = format_contact_id(get(dict, "doped_contact_id", -1))
    any(i -> i < 1, inactive_contact_id) && throw(ConfigFileError("Invalid doped_contact_id: missing or misspelled key"))
    bulk_imp_model = haskey(dict, "bulk_impurity_density") ? ImpurityDensity(T, dict["bulk_impurity_density"], input_units) : ConstantImpurityDensity{T}(-1e16)
    PtypePNJunctionImpurityDensity{T}(lithium_annealing_temperature, lithium_annealing_time, contact_with_lithium_doped, inactive_contact_id, bulk_imp_model)
end

function get_impurity_density(PtypePNjunction::PtypePNJunctionImpurityDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    get_impurity_density(PtypePNjunction.bulk_imp_model, pt)+get_impurity_density(PtypePNjunction.surface_imp_model, pt)
end

(*)(scale::Real, pidm::PtypePNJunctionImpurityDensity{T}) where {T} = PtypePNJunctionImpurityDensity{T}(scale * pidm.surface_imp_model, scale * pidm.bulk_imp_model)
(+)(offset::Union{<:Real, <:Quantity{<:Real, Unitful.𝐋^(-3)}}, pidm::PtypePNJunctionImpurityDensity{T}) where {T} = PtypePNJunctionImpurityDensity{T}(offset + pidm.surface_imp_model, offset + pidm.bulk_imp_model)

