abstract type AbstractSemiconductor{T} <: AbstractObject{T} end

struct Semiconductor{T,G,MT,CDM,IDM} <: AbstractSemiconductor{T}
    temperature::T
    material::MT
    impurity_density_model::IDM
    charge_drift_model::CDM
    geometry::G
end

function Semiconductor{T}(dict::Dict, input_units::NamedTuple, outer_transformations) where T <: SSDFloat
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
        cdm = getfield(Main, Symbol(dict["charge_drift_model"]["model"])){T}
        cdm <: AbstractChargeDriftModel{T} ? cdm() : ElectricFieldChargeDriftModel{T}()
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