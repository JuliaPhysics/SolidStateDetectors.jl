"""
struct ThermalDiffusionLithiumDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

Lithium impurity density model

ref: [Dai _et al._ (2023)](https://doi.org/10.1016/j.apradiso.2022.110638)
 
## Fields
* `lithium_annealing_temperature::T`: lithium annealing temperature in Kelvin, when the lithium is diffused into the crystal. The default value is 623 K.
* `lithium_annealing_time::T`: lithium annealing time in seconds. The default value is 18 minutes.
* `distance_to_contact::Function`: the function for describing the depth to surface. 
1) use the built-in function `ConstructiveSolidGeometry.distance_to_surface` to calculate the distance to the surface of the doped contact 
2) custom. 
Custom function might be much faster while the detector has good symmetry.
* `lithium_density_on_contact::T`: the lithium concentration in the surface (in m^-3)
* `lithium_diffusivity_in_germanium::T`: the diffusivity of lithium in germanium (in m^2*s^-1)
"""

struct ThermalDiffusionLithiumDensity{T <: SSDFloat, G <: Union{<:AbstractGeometry, Nothing}} <: AbstractImpurityDensity{T}
    lithium_annealing_temperature::T
    lithium_annealing_time::T
    distance_to_contact::Function
    lithium_density_on_contact::T
    lithium_diffusivity_in_germanium::T
end

function calculate_lithium_diffusivity_in_germanium(lithium_annealing_temperature::T)::T where {T <: SSDFloat}
    # D0 [m^2*s^-1]
    # H [cal]
    D0::T, H::T = ifelse(lithium_annealing_temperature <= 873, T.((2.5e-3*1e-4, 11800)), T.((1.3e-3*1e-4, 10700)))
    D0 * exp(-H/(R_gas*lithium_annealing_temperature))
end
function calculate_lithium_saturated_density(lithium_annealing_temperature::T)::T where {T <: SSDFloat}
    exp10(21.27 - 2610.0/(lithium_annealing_temperature) + 6)
end

function ThermalDiffusionLithiumDensity{T}(
    lithium_annealing_temperature::T,
    lithium_annealing_time::T,
    contact_with_lithium_doped::G;
    distance_to_contact::Function = pt::AbstractCoordinatePoint{T} -> ConstructiveSolidGeometry.distance_to_surface(pt, contact_with_lithium_doped),
    lithium_density_on_contact::T = calculate_lithium_saturated_density(lithium_annealing_temperature),
    lithium_diffusivity_in_germanium::T = calculate_lithium_diffusivity_in_germanium(lithium_annealing_temperature),
) where {T <: SSDFloat, G <: Union{<:AbstractGeometry, Nothing}}
    ThermalDiffusionLithiumDensity{T, G}(lithium_annealing_temperature, lithium_annealing_time, distance_to_contact, lithium_density_on_contact, lithium_diffusivity_in_germanium)
end

function ImpurityDensity(T::DataType, t::Val{:li_diffusion}, dict::AbstractDict, input_units::NamedTuple)
    lithium_annealing_temperature = _parse_value(T, get(dict, "lithium_annealing_temperature", 623u"K"), input_units.temperature)
    lithium_annealing_time = _parse_value(T, get(dict, "lithium_annealing_time", 18u"minute"), internal_time_unit)
    contact_with_lithium_doped = haskey(dict, "contact_with_lithium_doped") ? dict["contact_with_lithium_doped"] : nothing # you don't have to pass the geometry of doped contact only when the distance_to_contact is passed
    ThermalDiffusionLithiumDensity{T}(lithium_annealing_temperature, lithium_annealing_time, contact_with_lithium_doped)
end

function get_impurity_density(li_diffusion::ThermalDiffusionLithiumDensity{T, G}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat, G <: Union{<:AbstractGeometry, Nothing}}
    pt::CartesianPoint{T} = CartesianPoint(pt)
    depth = li_diffusion.distance_to_contact(pt)
    li_diffusion.lithium_density_on_contact * SpecialFunctions.erfc(depth/2/sqrt(li_diffusion.lithium_diffusivity_in_germanium*li_diffusion.lithium_annealing_time))
end
