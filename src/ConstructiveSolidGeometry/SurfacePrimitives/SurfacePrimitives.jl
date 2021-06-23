abstract type AbstractPlanarSurfacePrimitive{T} <: AbstractSurfacePrimitive{T} end
abstract type AbstractCurvedSurfacePrimitive{T} <: AbstractSurfacePrimitive{T} end

include("Plane.jl")
include("Polygon.jl")
include("ConeMantle.jl")
include("EllipticalSurface.jl")
include("EllipsoidMantle.jl")

