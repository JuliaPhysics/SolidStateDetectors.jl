# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

include("Units.jl")

include("Points.jl")
include("Vectors.jl")

abstract type AbstractGeometry{T, N} end
abstract type AbstractVolumePrimitive{T, N} <: AbstractGeometry{T, N} end
abstract type AbstractSurfacePrimitive{T, N} <: AbstractGeometry{T, N} end
abstract type AbstractSet{T, N} <: AbstractGeometry{T, N} end

function Geometry(T, geom_sub_dict::Dict, geom_unit)
    haskey(geom_sub_dict, "positive") ? geometry_positive = Geometries(T, geom_sub_dict["positive"], geom_unit) : geometry_positive = AbstractGeometry[]
    haskey(geom_sub_dict, "negative") ? geometry_negative = Geometries(T, geom_sub_dict["negative"], geom_unit) : geometry_negative = AbstractGeometry[]
    geometry = geometry_positive[1]
    for g in sort!(vcat(geometry_positive,geometry_negative))
        g in geometry_negative ? geometry -= g : nothing
        g in geometry_positive ? geometry += g : nothing
    end
    return geometry, geometry_positive, geometry_negative
end

function Geometry(T, geom_sub_dict::Vector, geom_unit)
    geometry_positive = Geometries(T, geom_sub_dict, geom_unit)
    geometry_negative = AbstractGeometry[]
    geometry = geometry_positive[1]
    for g in sort!(vcat(geometry_positive,geometry_negative))
        g in geometry_negative ? geometry -= g : nothing
        g in geometry_positive ? geometry += g : nothing
    end
    return geometry, geometry_positive, geometry_negative
end

function in(pt::AbstractCoordinatePoint{T}, vg::Vector{AbstractGeometry{T}})::Bool where {T <: SSDFloat}
    inside::Bool = false
    for g in vg
        if in(pt, g)
            inside = true
            break
        end
    end
    return inside
end

include("VolumePrimitives/VolumePrimitives.jl")
include("SurfacePrimitives/SurfacePrimitives.jl")
include("LinePrimitives/LinePrimitives.jl")

include("Sets.jl")


include("IO.jl")
