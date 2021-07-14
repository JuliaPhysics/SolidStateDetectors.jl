abstract type AbstractPassive{T} <: AbstractObject{T} end

mutable struct Passive{T,G,MT,CDM} <: AbstractPassive{T}
    name::String
    id::Int
    potential::T
    temperature::T
    material::MT
    charge_density_model::CDM
    geometry::G
end

function Passive{T}(dict::Dict, input_units::NamedTuple, outer_transformations) where T <: SSDFloat
    name = haskey(dict, "name") ? dict["name"] : "external part"
    id::Int = haskey(dict, "id") ? dict["id"] : -1
    potential = haskey(dict, "potential") ? T(dict["potential"]) : :floating
    material = material_properties[materials[dict["material"]]]
    temperature = haskey(dict, "temperature") ? T(dict["temperature"]) : T(293)
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
    println("\t________"*"Passive{$T} $(d.name) $(d.id)"*"________\n")
    println("\t---General Properties---")
    println("\t-Potential: \t\t $(d.potential) V")
    println("\t-Material: \t $(d.material.name)")
    println()
end
print(io::IO, d::Passive{T}) where {T} = print(io, "Passive $(d.name) - id $(d.id) - $(d.potential) V")
show(io::IO, d::Passive) = print(io, d)
show(io::IO,::MIME"text/plain", d::Passive) = show(io, d)
