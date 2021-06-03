abstract type AbstractPassive{T} <: AbstractObject{T} end

mutable struct Passive{T} <: AbstractPassive{T}
    name::String
    id::Int
    potential::Union{Symbol,T}
    temperature::Union{T,Missing}
    material::NamedTuple
    charge_density_model::AbstractChargeDensity{T}
    geometry::AbstractGeometry{T}
    geometry_positive::Vector{AbstractGeometry{T}}
    geometry_negative::Vector{AbstractGeometry{T}}
    decomposed_surfaces::Vector{AbstractGeometry{T}}
    
    Passive{T}() where T <: SSDFloat = new{T}()
end

function Passive{T}(dict::Dict, input_units::NamedTuple, transformations::Vector{CSGTransformation}) where T <: SSDFloat
    pass = Passive{T}()
    haskey(dict, "name") ? pass.name = dict["name"] : pass.name = "external part"
    haskey(dict, "id") ? pass.id = dict["id"] : pass.id = -1
    haskey(dict,"potential") ? pass.potential = T(dict["potential"]) : pass.potential = :floating
    haskey(dict, "temperature") ? pass.temperature = T(dict["temperature"]) : pass.temperature = missing
    pass.material = material_properties[materials[dict["material"]]]
    pass.charge_density_model = if haskey(dict, "charge_density") 
        haskey(dict, "charge_density") 
    elseif haskey(dict, "charge_density_model") 
        @warn "Config file deprication: There was an internal change from v0.5.1 to v0.6.0 regarding the 
            charge density of `Passive` objects. 
            Since v0.5.0, the elementary charge is not automatically multiplied to the distribution as it
            is a charge density and not an impurity density. The values in the config files should be adapted
            and the name of the field should be changed from \"charge_density_model\" into \"charge_density\".
            This warning will result in an error in later versions."
       haskey(dict, "charge_density_model") 
    else
        ConstantChargeDensity{T}(0)
    end
    pass.geometry = transform(Geometry(T, dict["geometry"], input_units), transformations)
    pass.geometry_positive, pass.geometry_negative = get_decomposed_volumes(pass.geometry)
    pass.decomposed_surfaces = vcat(get_decomposed_surfaces.(pass.geometry_positive)...)
    return pass
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
