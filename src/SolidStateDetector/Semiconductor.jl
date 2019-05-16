abstract type AbstractSemiconductor{T} <: AbstractObject{T} end

mutable struct Semiconductor{T} <: AbstractSemiconductor{T}
    name::String
    id::Int
    temperature::T
    material::NamedTuple
    bulk_type::Symbol
    charge_density_model::AbstractChargeDensityModel{T}
    geometry::AbstractGeometry{T}
    geometry_positive::Vector{AbstractGeometry{T}}
    geometry_negative::Vector{AbstractGeometry{T}}

    Semiconductor{T}() where T <: SSDFloat = new{T}()
end

function Semiconductor{T}(dict::Dict, inputunit_dict::Dict{String,Unitful.Units}) where T <: SSDFloat
    sc = Semiconductor{T}()
    sc.charge_density_model = if haskey(dict, "charge_density_model") 
        ChargeDensityModel(T, dict["charge_density_model"], inputunit_dict)
    else
        ZeroChargeDensityModel{T}()
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


function show(io::IO, d::Semiconductor{T}) where {T <: SSDFloat} println(d) end
function print(io::IO, d::Semiconductor{T}) where {T <: SSDFloat} println(d) end
function display(io::IO, d::Semiconductor{T} ) where {T <: SSDFloat} println(d) end
function show(io::IO,::MIME"text/plain", d::Semiconductor) where {T <: SSDFloat}
    show(io, d)
end
