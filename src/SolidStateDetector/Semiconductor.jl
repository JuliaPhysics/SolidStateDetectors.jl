abstract type AbstractSemiconductor{T} <: AbstractObject{T} end

struct Semiconductor{T,G,MT,CDM,IDM} <: AbstractSemiconductor{T}
    temperature::T
    material::MT
    impurity_density_model::IDM
    charge_drift_model::CDM
    geometry::G
end

function Semiconductor{T}(dict::Dict, input_units::NamedTuple, transformations = missing) where T <: SSDFloat
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
    temperature = haskey(dict, "temperature") ? T(dict["temperature"]) : T(80)
    material = material_properties[materials[dict["material"]]]
    geometry = transform(Geometry(T, dict["geometry"], input_units), transformations)
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
