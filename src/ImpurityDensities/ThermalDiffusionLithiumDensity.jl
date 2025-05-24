"""
    struct ThermalDiffusionLithiumDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

Lithium impurity density model

ref: [Dai _et al._ (2023)](https://doi.org/10.1016/j.apradiso.2022.110638)
 
## Fields
* `Li_annealing_temperature::T`: lithium annealing temperature when the lithium is diffused into the crystal. The default value is 623 K.
* `Li_annealing_time::T`: lithium annealing time. The default value is 18 minutes.
* `calculate_depth2surface::Function`: the function for describing the depth to surface
* `Li_surface::T`: the lithium concentration in the surface
* `D_Li::T`: the diffusivity of lithium in germanium
"""

struct ThermalDiffusionLithiumDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}
	Li_annealing_temperature::T # in K
	Li_annealing_time::T # in s
    calculate_depth2surface::Function # pt -> depth (depth in m)
	Li_surface::T # in m^-3
    D_Li::T # in m^2*s^-1
end
function calculate_D_Li(Li_annealing_temperature::T)::T where {T <: SSDFloat}
    if Li_annealing_temperature<=873
        D0 = 2.5e-3*1e-4 # m^2*s^-1
        H = 11800 # cal
    else
        D0 = 1.3e-3*1e-4
        H = 10700
    end
    D_Li = D0 * exp(-H/(R_gas*Li_annealing_temperature))
    D_Li
end
function calculate_Li_surface_saturated_density(Li_annealing_temperature::T)::T where {T <: SSDFloat}
    10^(21.27 - 2610.0/(Li_annealing_temperature)) * 1e6 # m^-3
end
function ThermalDiffusionLithiumDensity{T}(Li_annealing_temperature::T,
    Li_annealing_time::T,
    calculate_depth2surface::Function,
    Li_surface::T=calculate_Li_surface_saturated_density(Li_annealing_temperature),
    D_Li::T=calculate_D_Li(Li_annealing_temperature)
    ) where {T <: SSDFloat}
    ThermalDiffusionLithiumDensity(Li_annealing_temperature,Li_annealing_time,calculate_depth2surface, Li_surface, D_Li)
end

function ImpurityDensity(T::DataType, t::Val{:tdld}, dict::AbstractDict, input_units::NamedTuple)
    Li_annealing_temperature = T(623)
    Li_annealing_time= T(18*60)
    calculate_depth2surface = pt->0.001

    if haskey(dict, "Li_annealing_temperature")
        Li_annealing_temperature = _parse_value(T, dict["Li_annealing_temperature"], u"K")
    end
    if haskey(dict, "Li_annealing_time")
        Li_annealing_time = _parse_value(T, dict["Li_annealing_time"], u"s")
    end
    if haskey(dict, "calculate_depth2surface")
        calculate_depth2surface = dict["calculate_depth2surface"]
    end

    ThermalDiffusionLithiumDensity{T}(Li_annealing_temperature, Li_annealing_time, calculate_depth2surface)
end

function get_impurity_density(tdld::ThermalDiffusionLithiumDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    pt::CartesianPoint{T} = CartesianPoint(pt)
    depth = tdld.calculate_depth2surface(pt)
	tdld.Li_surface * Distributions.erfc(depth/2/sqrt(tdld.D_Li*tdld.Li_annealing_time))
end