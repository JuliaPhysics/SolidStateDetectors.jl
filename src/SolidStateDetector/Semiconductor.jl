abstract type AbstractSemiconductor{T} <: AbstractObject{T} end

"""
    struct Semiconductor{T,G,MT,CDM,IDM} <: AbstractSemiconductor{T}
        
Semiconductor bulk of a [`SolidStateDetector`](@ref).

This is the volume in which electrons and holes will drift during the signal development.

## Parametric types
* `T`: Precision type.
* `G`: Type of `geometry`.
* `MT`: Type of `material`.
* `CDM`: Type of `charge_drift_model`.
* `IDM`: Type of `impurity_density_model`.

## Fields
* `temperature::T`: Temperature (in K) of the semiconductor.
* `material::MT`: Material of the semiconductor.
* `impurity_density_model::IDM`: Impurity density model for the points inside the semiconductor.
* `charge_drift_model::CDM`: Model that describes the drift of electrons and holes inside the semiconductor.
* `geometry::G`: Geometry of the semiconductor, see [Constructive Solid Geometry (CSG)](@ref).

## Definition in Configuration File

A `Semiconductor` is defined through the `semiconductor` field of a detector.
It needs `material` and `geometry`, and can optionally be given a `temperature`, `impurity_density` and `charge_drift_model`.

An example definition of a semiconductor looks like this:
```yaml 
semiconductor:
  material: HPGe
  temperature: 78
  impurity_density: # ...
  charge_drift_model: # ...
  geometry: # ...
```
"""
struct Semiconductor{T,G,MT,CDM,IDM} <: AbstractSemiconductor{T}
    temperature::T
    material::MT
    impurity_density_model::IDM
    charge_drift_model::CDM
    geometry::G
end

function Semiconductor{T}(dict::AbstractDict, input_units::NamedTuple, outer_transformations) where T <: SSDFloat
    impurity_density_model = if haskey(dict, "impurity_density") 
        ImpurityDensity(T, dict["impurity_density"], input_units)
    elseif haskey(dict, "charge_density_model") 
        @warn "Config file deprication: The field \"charge_density_model\" under semiconductor is deprecated. 
            It should be changed to \"impurity_density\". In later version this will result in an error.
            For now, it will be treated as an impurity density."
        ImpurityDensity(T, dict["charge_density_model"], input_units)
    else
        ConstantImpurityDensity{T}(0)
    end
    charge_drift_model = if haskey(dict, "charge_drift_model") && haskey(dict["charge_drift_model"], "model")
        model = Symbol(dict["charge_drift_model"]["model"])
        cdm = if model in names(SolidStateDetectors, all = true) && getfield(SolidStateDetectors, model) <: AbstractChargeDriftModel
            getfield(SolidStateDetectors, model){T}()
        else
            throw(ConfigFileError("There is no charge drift model called `$(dict["charge_drift_model"]["model"])`."))
        end
    else
        ElectricFieldChargeDriftModel{T}()
    end
    material = material_properties[materials[dict["material"]]]
    temperature = if haskey(dict, "temperature") 
        T(dict["temperature"])
    elseif material.name == "High Purity Germanium"
        T(78)
    else
        T(293)
    end
    inner_transformations = parse_CSG_transformation(T, dict, input_units)
    transformations = combine_transformations(inner_transformations, outer_transformations)
    geometry = Geometry(T, dict["geometry"], input_units, transformations)
    return Semiconductor(temperature, material, impurity_density_model, charge_drift_model, geometry)
end

function println(io::IO, d::Semiconductor{T}) where {T <: SSDFloat}
    println("\t---General Properties---")
    println("\t-Detector Material: \t $(d.material.name)")
    println("\t-Charge Drift Model:\t $(typeof(d.charge_drift_model).name.name)")
end

print(io::IO, d::Semiconductor{T}) where {T} = print(io, "Semiconductor{$T} - $(d.material.name)")

show(io::IO, d::Semiconductor) = print(io, d)
show(io::IO,::MIME"text/plain", d::Semiconductor) = show(io, d)


function Semiconductor(sc::Semiconductor{T,G,MT,CDM,IDM}, impurity_density::AbstractImpurityDensity{T}) where {T,G,MT,CDM,IDM}
    Semiconductor(sc.temperature, sc.material, impurity_density, sc.charge_drift_model, sc.geometry)
end
function Semiconductor(sc::Semiconductor{T,G,MT,CDM,IDM}, chargedriftmodel::AbstractChargeDriftModel{T}) where {T,G,MT,CDM,IDM}
    Semiconductor(sc.temperature, sc.material, sc.impurity_density_model, chargedriftmodel, sc.geometry)
end

"""
    scaling_factor_for_permittivity_in_undepleted_region(sc::Semiconductor{T})::T where {T}

This function is called in the calculations of weighting potentials of undepleted detectors. 
The electric permittivity, ``ϵ_{r}``, is scaled with this function in areas where the detector is undepleted.
A value between `[0, +Inf]` should be returned. However, `Inf` should not be returned but instead
a very high value should be returned in order to mimic perfect conductivity if that is desired. 

## Arguments
* `sc::Semiconductor{T}`: Semiconductor for which the dielectric permittivity should be scaled up.

!!! danger "Experimental feature!"
    This feature is under research. The goal is to study the properties / signal response of undepleted detector. 
    This function is indented to be overwritten by the user to study the response. 
"""
function scaling_factor_for_permittivity_in_undepleted_region(sc::Semiconductor{T})::T where {T}
    100000
end

