abstract type AbstractFlatSurfacePrimitive{T} <: AbstractSurfacePrimitive{T} end
abstract type AbstractBentSurfacePrimitive{T} <: AbstractSurfacePrimitive{T} end

include("Plane.jl")
include("Polygon.jl")
include("ConeMantle.jl")
include("EllipticalSurface.jl")