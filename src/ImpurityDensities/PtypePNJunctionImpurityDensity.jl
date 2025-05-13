"""
    struct ThermalDiffusionLithiumDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

A PN junctin impurity model based on lithium thermal diffusion and constant bulk impurity density.
The surface lithium density is saturated.

ref: [Dai _et al._ (2023)](https://doi.org/10.1016/j.apradiso.2022.110638)
 
## Fields
* `Li_annealing_temperature::T`: lithium annealing temperature when the lithium is diffused into the crystal. The default value is 623 K.
* `Li_annealing_time::T`: lithium annealing time. The default value is 18 minutes.
* `calculate_depth2surface::Function`: the function for describing the depth to surface.
* `bulk_imp::T`: the bulk impurity density constant.
* `bulk/surface_imp_model::T`: there two variables will be constructed with above parameters.
"""

struct PtypePNJunctionImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}
	Li_annealing_temperature::T
	Li_annealing_time::T
    calculate_depth2surface::Function
    bulk_imp::T
    surface_imp_model::ThermalDiffusionLithiumDensity{T}
    bulk_imp_model::ConstantImpurityDensity{T}
end
function PtypePNJunctionImpurityDensity{T}(Li_annealing_temperature::T,
    Li_annealing_time::T,
    calculate_depth2surface::Function,
    bulk_imp::T
    ) where {T <: SSDFloat}
    bulk_imp_model = ConstantImpurityDensity{T}(bulk_imp)
    surface_imp_model = ThermalDiffusionLithiumDensity{T}(Li_annealing_temperature, Li_annealing_time, calculate_depth2surface)
    PtypePNJunctionImpurityDensity(Li_annealing_temperature,Li_annealing_time,calculate_depth2surface, bulk_imp, surface_imp_model, bulk_imp_model)
end

function ImpurityDensity(T::DataType, t::Val{:PtypePNjunction}, dict::AbstractDict, input_units::NamedTuple)
    Li_annealing_temperature = T(623)
    Li_annealing_time= T(18*60)
    calculate_depth2surface = pt->0.01-hypot(pt[1],pt[2]) # default for the public_TrueCoaxial.yaml
    bulk_imp = T(-1e16)

    if haskey(dict, "Li_annealing_temperature")
        Li_annealing_temperature = _parse_value(T, dict["Li_annealing_temperature"], u"K")
    end
    if haskey(dict, "Li_annealing_time")
        Li_annealing_time = _parse_value(T, dict["Li_annealing_time"], u"s")
    end
    if haskey(dict, "calculate_depth2surface")
        calculate_depth2surface = dict["calculate_depth2surface"]
    end
    if haskey(dict, "bulk_imp")
        bulk_imp = _parse_value(T, dict["bulk_imp"], input_units.length^(-3))
    end

    PtypePNJunctionImpurityDensity{T}(Li_annealing_temperature, Li_annealing_time, calculate_depth2surface, bulk_imp)
end

function get_impurity_density(PtypePNjunction::PtypePNJunctionImpurityDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    get_impurity_density(PtypePNjunction.bulk_imp_model, pt)+get_impurity_density(PtypePNjunction.surface_imp_model, pt)
end