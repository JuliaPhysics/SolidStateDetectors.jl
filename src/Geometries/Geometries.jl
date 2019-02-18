# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

include("Units.jl")

include("Points.jl")
include("Vectors.jl")

abstract type AbstractGeometry{T, N} end
abstract type AbstractVolumePrimitive{T, N} <: AbstractGeometry{T, N} end
abstract type AbstractSurfacePrimitive{T, N} <: AbstractGeometry{T, N} end
abstract type AbstractSet{T, N} <: AbstractGeometry{T, N} end

include("LinePrimitives/LinePrimitives.jl")
include("SurfacePrimitives/SurfacePrimitives.jl")
include("VolumePrimitives/VolumePrimitives.jl")

include("Sets.jl")


include("IO.jl")