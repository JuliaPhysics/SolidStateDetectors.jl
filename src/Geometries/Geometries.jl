# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

include("Units.jl")

include("Points.jl")
include("Vectors.jl")

abstract type AbstractGeometry{T, N} end
abstract type AbstractVolumePrimitive{T, N} <: AbstractGeometry{T, N} end
abstract type AbstractSurfacePrimitive{T, N} <: AbstractGeometry{T, N} end
abstract type AbstractSet{T, N} <: AbstractGeometry{T, N} end


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
