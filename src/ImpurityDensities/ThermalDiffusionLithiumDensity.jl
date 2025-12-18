"""
struct ThermalDiffusionLithiumDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

Lithium impurity density model. Ref: [Dai _et al._ (2023)](https://doi.org/10.1016/j.apradiso.2022.110638)
 
## Fields
* `lithium_annealing_temperature::T`: lithium annealing temperature in Kelvin, when the lithium is diffused into the crystal. The default value is 623 K.
* `lithium_annealing_time::T`: lithium annealing time in seconds. The default value is 18 minutes.
* `distance_to_contact::Function`: the function for describing the depth to surface. 
    1) use the built-in function `ConstructiveSolidGeometry.distance_to_surface` to calculate the distance to the surface of the doped contact 
    2) custom. Custom function might be much faster while the detector has good symmetry.
* `lithium_density_on_contact::T`: the lithium concentration in the surface (in m^-3)
* `lithium_diffusivity_in_germanium::T`: the diffusivity of lithium in germanium (in m^2*s^-1)

## Parameters for constructing the model
* `lithium_annealing_temperature::T`
* `lithium_annealing_time::T`
* `contact_with_lithium_doped::G`: the geometry of the contact with lithium doped, which is used to calculate the distance to the surface.
    - If you don't pass this parameter, the `distance_to_contact` function should be passed.
    - If you pass this parameter but don't pass the `distance_to_contact` function, the `distance_to_contact` function will use the default one.
* `distance_to_contact::Function`: optional, default is `ConstructiveSolidGeometry.distance_to_surface(pt, contact_with_lithium_doped)`
* `lithium_density_on_contact::T`: optional, default is the saturated lithium density at the given temperature
* `lithium_diffusivity_in_germanium::T`: optional, default is calculated with the annealing temperature
"""
struct ThermalDiffusionLithiumDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}
    lithium_annealing_temperature::T
    lithium_annealing_time::T
    inactive_contact_id::Int
    distance_to_contact::Function
    lithium_density_on_contact::T
    lithium_diffusivity_in_germanium::T
end

struct LithiumDiffusionBranch{T <: SSDFloat}
    T_min::T
    T_max::T
    D0::T
    H::T
end

struct LithiumDiffusionParameters{T <: SSDFloat}
    lowT::LithiumDiffusionBranch{T}
    highT::LithiumDiffusionBranch{T}
end

struct LithiumSaturationParameters{T <: SSDFloat}
    a::T
    b::T
end

struct ThermalDiffusionLithiumDensityParameters{T <: SSDFloat}
    diffusion::LithiumDiffusionParameters{T}
    saturation::LithiumSaturationParameters{T}
end

include("ThermalDiffusionLithiumDensityParameters.jl")

function calculate_lithium_diffusivity_in_germanium(lithium_annealing_temperature::T,parameters::LithiumDiffusionParameters)::T where {T <: SSDFloat}
    # D0 [m^2*s^-1]
    # H [cal]
    D0::T, H::T = ifelse(lithium_annealing_temperature <= parameters.lowT.T_max, T.((parameters.lowT.D0, parameters.lowT.H)), T.((parameters.highT.D0, parameters.highT.H)))
    D0 * exp(-H/(R_gas*lithium_annealing_temperature))
end
function calculate_lithium_saturated_density(lithium_annealing_temperature::T,parameters::LithiumSaturationParameters)::T where {T <: SSDFloat}
    exp10(parameters.a - parameters.b/lithium_annealing_temperature)
end

function ThermalDiffusionLithiumDensity{T}( 
    lithium_annealing_temperature::T,
    lithium_annealing_time::T,
    contact_with_lithium_doped::G,
    inactive_contact_id::Int;
    modelParameters::ThermalDiffusionLithiumDensityParameters{T} = ThermalDiffusionLithiumParameters(), 
    distance_to_contact::Function = pt::AbstractCoordinatePoint{T} -> ConstructiveSolidGeometry.distance_to_surface(pt, contact_with_lithium_doped),
    lithium_density_on_contact::T = calculate_lithium_saturated_density(lithium_annealing_temperature,modelParameters.saturation),
    lithium_diffusivity_in_germanium::T = calculate_lithium_diffusivity_in_germanium(lithium_annealing_temperature,modelParameters.diffusion),
) where {T <: SSDFloat, G <: Union{<:AbstractGeometry, Nothing}}
    ThermalDiffusionLithiumDensity{T}(lithium_annealing_temperature, lithium_annealing_time, inactive_contact_id, distance_to_contact, lithium_density_on_contact, lithium_diffusivity_in_germanium)
end

function ImpurityDensity(T::DataType, t::Val{:li_diffusion}, dict::AbstractDict, input_units::NamedTuple)
    lithium_annealing_temperature = _parse_value(T, get(dict, "lithium_annealing_temperature", 623u"K"), input_units.temperature)
    lithium_annealing_time = _parse_value(T, get(dict, "lithium_annealing_time", 18u"minute"), internal_time_unit)
    contact_with_lithium_doped = haskey(dict, "contact_with_lithium_doped") ? dict["contact_with_lithium_doped"] : nothing # you don't have to pass the geometry of doped contact only when the distance_to_contact is passed
    inactive_contact_id = get(dict, "doped_contact_id", -1)
    inactive_contact_id < 1 && error("Invalid doped_contact_id: missing or misspelled key")
    ThermalDiffusionLithiumDensity{T}(lithium_annealing_temperature, lithium_annealing_time, contact_with_lithium_doped, inactive_contact_id)
end

function get_impurity_density(li_diffusion::ThermalDiffusionLithiumDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    pt::CartesianPoint{T} = CartesianPoint(pt)
    depth = li_diffusion.distance_to_contact(pt)
    li_diffusion.lithium_density_on_contact * SpecialFunctions.erfc(depth/2/sqrt(li_diffusion.lithium_diffusivity_in_germanium*li_diffusion.lithium_annealing_time))
end

(*)(scale::Real, tidm::ThermalDiffusionLithiumDensity{T}) where {T} = ThermalDiffusionLithiumDensity{T}(tidm.lithium_annealing_temperature, tidm.lithium_annealing_time, tidm.inactive_contact_id, tidm.distance_to_contact, T(scale * tidm.lithium_density_on_contact), tidm.lithium_diffusivity_in_germanium)

(+)(offset::Union{<:Real, <:Quantity{<:Real, Unitful.𝐋^(-3)}}, tidm::ThermalDiffusionLithiumDensity{T}) where {T} = ThermalDiffusionLithiumDensity{T}(tidm.lithium_annealing_temperature, tidm.lithium_annealing_time, tidm.inactive_contact_id, tidm.distance_to_contact, T(to_internal_units(offset) + tidm.lithium_density_on_contact), tidm.lithium_diffusivity_in_germanium)
