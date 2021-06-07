abstract type AbstractSemiconductor{T} <: AbstractObject{T} end

mutable struct Semiconductor{T} <: AbstractSemiconductor{T}
    name::String
    id::Int
    temperature::T
    material::NamedTuple
    impurity_density_model::AbstractImpurityDensity{T}
    charge_drift_model::AbstractChargeDriftModel{T}
    geometry::AbstractGeometry{T}
    geometry_positive::Vector{AbstractGeometry{T}}
    geometry_negative::Vector{AbstractGeometry{T}}
    decomposed_surfaces::Vector{AbstractGeometry{T}}

    Semiconductor{T}() where T <: SSDFloat = new{T}()
end

function Semiconductor{T}(dict::Dict, input_units::NamedTuple, transformations::Vector{CSGTransformation} = []) where T <: SSDFloat
    sc = Semiconductor{T}()
    sc.impurity_density_model = if haskey(dict, "impurity_density") 
        ImpurityDensity(T, dict["impurity_density"], input_units)
    elseif haskey(dict, "charge_density_model") 
        @warn "Config file deprication: The field \"charge_density_model\" under semiconductor is deprecated. 
            It should be changed to \"impurity_density\". In later version this will result in an error.
            For now, it will be treated as an impurity density."
        ImpurityDensity(T, dict["charge_density_model"], input_units)
    else
        ConstantImpurityDensity{T}(0)
    end
    sc.charge_drift_model = if haskey(dict, "charge_drift_model") && haskey(dict["charge_drift_model"], "model")
        cdm = getfield(Main, Symbol(dict["charge_drift_model"]["model"])){T}
        cdm <: AbstractChargeDriftModel{T} ? cdm() : ElectricFieldChargeDriftModel{T}()
    else
        ElectricFieldChargeDriftModel{T}()
    end
    sc.material = material_properties[materials[dict["material"]]]
    sc.geometry = transform(Geometry(T, dict["geometry"], input_units), transformations)
    sc.geometry_positive, sc.geometry_negative = get_decomposed_volumes(sc.geometry)
    sc.decomposed_surfaces = vcat(get_decomposed_surfaces.(sc.geometry_positive)...)
    return sc
end


function println(io::IO, d::Semiconductor{T}) where {T <: SSDFloat}
    println("\t---General Properties---")
    println("\t-Detector Material: \t $(d.material.name)")
    println("\t-Charge Drift Model:\t $(typeof(d.charge_drift_model).name.name)")
end

print(io::IO, d::Semiconductor{T}) where {T} = print(io, "Semiconductor{$T} - $(d.material.name)")

show(io::IO, d::Semiconductor) = print(io, d)
show(io::IO,::MIME"text/plain", d::Semiconductor) = show(io, d)
