abstract type AbstractPassive{T} <: AbstractObject{T} end

"""
    struct Passive{T,G,MT,CDM} <: AbstractPassive{T}
        
Passive object, assigned to a [`SolidStateDetector`](@ref).

For the calculation of the [`ElectricPotential`](@ref) and [`WeightingPotential`](@ref), 
passives can be fixed to a constant potential. They can additionally have a charge density 
profile that has an influence on the [`ElectricPotential`](@ref).

## Parametric types
* `T`: Precision type.
* `G`: Type of `geometry`.
* `MT`: Type of `material`.
* `CDM`: Type of `charge_density_model`.

## Fields
* `name::String`: Custom name for the passive, relevant for plotting.
* `id::Int`: Unique id that will unambiguously identify the passive.
* `potential::T`: Potential (in V) to which the passive will be fixed during the calculation of the electric potential.
    For floating passives, the `potential` value is `NaN`.
* `temperature::T`: Temperature (in K) of the passive.
* `material::MT`: Material of the passive.
* `charge_density_model::CDM`: Charge density model for the points inside the passive.
* `geometry::G`: Geometry of the passive, see [Constructive Solid Geometry (CSG)](@ref).

## Definition in Configuration File

A `Passive` is defined through an entry in the `passives` array of a detector
or an entry in the `surroundings` array in the configuration file.
It needs `material` and `geometry` and can optionally be given a `name`, `id`, `potential`, `temperature` and `charge_density`.

An example definition of passives looks like this:
```yaml 
passives:
  - name: Passivated Surface
    material: HPGe
    charge_density: # ...
    geometry: # ...
  - name: Cryostat
    id: 3
    potential: 0
    temperature: 293K
    material: Al
    geometry: # ...
```
"""
struct Passive{T,G,MT,CDM} <: AbstractPassive{T}
    name::String
    id::Int
    potential::T # NaN is floating
    temperature::T
    material::MT
    charge_density_model::CDM
    geometry::G
end

const POTENTIAL_FLOATING = NaN

function Passive{T}(dict::AbstractDict, input_units::NamedTuple, outer_transformations) where {T <: SSDFloat}
    name::String = get(dict, "name", "External part")
    id::Int = get(dict, "id", -1)
    potential::T = _parse_value(T, get(dict, "potential", POTENTIAL_FLOATING), input_units.potential)
    material = material_properties[materials[dict["material"]]]
    temperature::T = _parse_value(T, get(dict, "temperature", 293u"K"), input_units.temperature)
    charge_density_model = if haskey(dict, "charge_density") 
        ChargeDensity(T, dict["charge_density"], input_units)
    elseif haskey(dict, "charge_density_model") 
        @warn "Configuration file deprecation: There was an internal change from v0.5.3 to v0.6.0 regarding the 
            charge density of `Passive` objects. 
            Since v0.6.0, the elementary charge is not automatically multiplied to the distribution as it
            is a charge density and not an impurity density. The values in the config files should be adapted
            and the name of the field should be changed from \"charge_density_model\" into \"charge_density\".
            This warning will result in an error in later versions."
        ChargeDensity(T, dict["charge_density_model"], input_units)
    else
        ConstantChargeDensity{T}(0)
    end
    inner_transformations = parse_CSG_transformation(T, dict, input_units)
    transformations = combine_transformations(inner_transformations, outer_transformations)
    geometry = Geometry(T, dict["geometry"], input_units, transformations)
    return Passive(name, id, potential, temperature, material, charge_density_model, geometry)
end

function println(io::IO, d::Passive{T}) where {T}
    println("\t________"*"$(d.name)"*"________\n")
    println("\t---General Properties---")
    println("\t-Potential: \t\t " * (isnan(d.potential) ? "floating" :  "$(d.potential) V"))
    println("\t-Material:  \t\t $(d.material.name)")
    println()
end
print(io::IO, d::Passive{T}) where {T} = print(io, "Passive $(d.name) - id $(d.id) - $(d.potential) V")
show(io::IO, d::Passive) = print(io, d)
show(io::IO,::MIME"text/plain", d::Passive) = show(io, d)
