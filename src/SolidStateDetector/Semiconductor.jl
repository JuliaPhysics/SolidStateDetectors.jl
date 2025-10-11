abstract type AbstractSemiconductor{T} <: AbstractObject{T} end

"""
    struct Semiconductor{T,G,MT,CDM,IDM,CTM} <: AbstractSemiconductor{T}
        
Semiconductor bulk of a [`SolidStateDetector`](@ref).

This is the volume in which electrons and holes will drift during the signal development.

## Parametric types
* `T`: Precision type.
* `G`: Type of `geometry`.
* `MT`: Type of `material`.
* `CDM`: Type of `charge_drift_model`.
* `IDM`: Type of `impurity_density_model`.
* `CTM`: Type of `charge_trapping_model`.

## Fields
* `temperature::T`: Temperature (in K) of the semiconductor.
* `material::MT`: Material of the semiconductor.
* `impurity_density_model::IDM`: Impurity density model for the points inside the semiconductor.
* `charge_drift_model::CDM`: Model that describes the drift of electrons and holes inside the semiconductor.
* `charge_trapping_model::CTM`: Model that describes the trapping of electrons and holes inside the semiconductor.
* `geometry::G`: Geometry of the semiconductor, see [Constructive Solid Geometry (CSG)](@ref).

## Definition in Configuration File

A `Semiconductor` is defined through the `semiconductor` field of a detector.
It needs `material` and `geometry`, and can optionally be given a `temperature`, `impurity_density`, `charge_drift_model` and `charge_trapping_model`.

An example definition of a semiconductor looks like this:
```yaml 
semiconductor:
  material: HPGe
  temperature: 78
  impurity_density: # ...
  charge_drift_model: # ...
  charge_trapping_model: # ...
  geometry: # ...
```
"""
struct Semiconductor{T,G,MT,CDM,IDM,CTM} <: AbstractSemiconductor{T}
    temperature::T
    material::MT
    impurity_density_model::IDM
    charge_drift_model::CDM
    charge_trapping_model::CTM
    geometry::G
end

function Semiconductor{T}(dict::AbstractDict, input_units::NamedTuple, outer_transformations) where T <: SSDFloat

    impurity_density_model = if haskey(dict, "impurity_density") 
        hascorrections = haskey(dict["impurity_density"], "corrections")
        impurity_density_scale = hascorrections ? _parse_value(T, get(dict["impurity_density"]["corrections"], "scale", 1), NoUnits) : T(1)
        impurity_density_offset = hascorrections ? _parse_value(T, get(dict["impurity_density"]["corrections"], "offset", 0), input_units.length^(-3)) : T(0)
        impurity_density_scale * ImpurityDensity(T, dict["impurity_density"], input_units) + impurity_density_offset
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
        cdm = if isdefined(SolidStateDetectors, model) && getfield(SolidStateDetectors, model) <: AbstractChargeDriftModel
            if model == :InactiveLayerChargeDriftModel
                InactiveLayerChargeDriftModel{T}(dict["charge_drift_model"], impurity_density_model, input_units)
            else
                getfield(SolidStateDetectors, model){T}(dict["charge_drift_model"])
            end
        else
            throw(ConfigFileError("There is no charge drift model called `$(dict["charge_drift_model"]["model"])`."))
        end
    else
        ElectricFieldChargeDriftModel{T}()
    end

    material = material_properties[materials[dict["material"]]]
    temperature = if haskey(dict, "temperature") 
        _parse_value(T, dict["temperature"], input_units.temperature)
    elseif material.name == "High Purity Germanium"
        T(78)
    else
        T(293)
    end

    inner_transformations = parse_CSG_transformation(T, dict, input_units)
    transformations = combine_transformations(inner_transformations, outer_transformations)
    geometry = Geometry(T, dict["geometry"], input_units, transformations)


    ctm_dict = haskey(dict, "charge_trapping_model") ? deepcopy(dict["charge_trapping_model"]) : Dict{String, Any}()

    if haskey(ctm_dict, "inactive_layer_geometry")
        ctm_dict["inactive_layer_geometry"] = Geometry(T, ctm_dict["inactive_layer_geometry"], input_units, transformations)
    end
    
    charge_trapping_model = if haskey(ctm_dict, "model_inactive")
        if haskey(ctm_dict, "parameters_inactive") && haskey(ctm_dict, "parameters") &&
            ctm_dict["parameters"]==ctm_dict["parameters_inactive"] && ctm_dict["model"] == "ConstantLifetime"
            ConstantLifetimeChargeTrappingModel{T}(ctm_dict)
        else
            CombinedChargeTrappingModel{T}(ctm_dict, temperature = temperature)
        end
        
    elseif haskey(ctm_dict, "model") && !haskey(ctm_dict, "model_inactive") && ctm_dict["model"] == "Boggs"
        BoggsChargeTrappingModel{T}(ctm_dict, temperature = temperature)
        
    elseif haskey(ctm_dict, "model") && !haskey(ctm_dict, "model_inactive") && ctm_dict["model"] == "ConstantLifetime"
        ConstantLifetimeChargeTrappingModel{T}(ctm_dict)

    else
        NoChargeTrappingModel{T}()
    end
    
    return Semiconductor(temperature, material, impurity_density_model, charge_drift_model, charge_trapping_model, geometry)
end

function println(io::IO, d::Semiconductor{T}) where {T <: SSDFloat}
    println("\t---General Properties---")
    println("\t-Detector Material: \t $(d.material.name)")
    println("\t-Impurity Density Model: $(typeof(d.impurity_density_model).name.name)")
    println("\t-Charge Drift Model:\t $(typeof(d.charge_drift_model).name.name)")
    if !(d.charge_trapping_model isa NoChargeTrappingModel)
        println("\t-Charge Trapping Model:\t $(typeof(d.charge_trapping_model).name.name)")
    end
end

print(io::IO, d::Semiconductor{T}) where {T} = print(io, "Semiconductor{$T} - $(d.material.name)")

show(io::IO, d::Semiconductor) = print(io, d)
show(io::IO,::MIME"text/plain", d::Semiconductor) = show(io, d)


function Semiconductor(sc::Semiconductor{T}, impurity_density::AbstractImpurityDensity{T}) where {T <: SSDFloat}
    Semiconductor(sc.temperature, sc.material, impurity_density, sc.charge_drift_model, sc.charge_trapping_model, sc.geometry)
end
function Semiconductor(sc::Semiconductor{T}, charge_drift_model::AbstractChargeDriftModel{T}) where {T <: SSDFloat}
    Semiconductor(sc.temperature, sc.material, sc.impurity_density_model, charge_drift_model, sc.charge_trapping_model, sc.geometry)
end
function Semiconductor(sc::Semiconductor{T}, charge_trapping_model::AbstractChargeTrappingModel{T}) where {T <: SSDFloat}
    Semiconductor(sc.temperature, sc.material, sc.impurity_density_model, sc.charge_drift_model, charge_trapping_model, sc.geometry)
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

