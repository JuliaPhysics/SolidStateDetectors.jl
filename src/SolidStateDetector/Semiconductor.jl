abstract type AbstractSemiconductor{T} <: AbstractObject{T} end

mutable struct Semiconductor{T} <: AbstractSemiconductor{T}
    name::String
    id::Int
    temperature::T
    material::NamedTuple
    bulk_type::Symbol
    charge_density_model::AbstractChargeDensity{T}
    geometry::AbstractGeometry{T}
    geometry_positive::Vector{AbstractGeometry{T}}
    geometry_negative::Vector{AbstractGeometry{T}}

    Semiconductor{T}() where T <: SSDFloat = new{T}()
end

function Semiconductor{T}(dict::Dict, inputunit_dict::Dict{String,Unitful.Units}) where T <: SSDFloat
    sc = Semiconductor{T}()
    sc.charge_density_model = if haskey(dict, "charge_density_model") 
        ChargeDensity(T, dict["charge_density_model"], inputunit_dict)
    else
        ZeroChargeDensity{T}()
    end
    sc.material = material_properties[materials[dict["material"]]]
    sc.bulk_type = bulk_types[dict["bulk_type"]]
    sc.geometry = Geometry(T, dict["geometry"], inputunit_dict)
    sc.geometry_positive, sc.geometry_negative = get_decomposed_volumes(sc.geometry)
    return sc
end


function println(io::IO, d::Semiconductor{T}) where {T <: SSDFloat}
    println("\t---General Properties---")
    println("\t-Detector Material: \t $(d.material.name)")
    println("\t-Bulk type: \t\t $(d.bulk_type)")
end

print(io::IO, d::Semiconductor{T}) where {T} = print(io, "Semiconductor{$T} - $(d.bulk_type) - $(d.material.name)")

show(io::IO, d::Semiconductor) = print(io, d)
show(io::IO,::MIME"text/plain", d::Semiconductor) = show(io, d)
