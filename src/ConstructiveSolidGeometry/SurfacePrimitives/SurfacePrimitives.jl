abstract type AbstractPlanarSurfacePrimitive{T} <: AbstractSurfacePrimitive{T} end
abstract type AbstractCurvedSurfacePrimitive{T} <: AbstractSurfacePrimitive{T} end

include("Plane.jl")
include("ConeMantle.jl")
include("EllipticalSurface.jl")
include("EllipsoidMantle.jl")
include("Polygon.jl")
include("TorusMantle.jl")

# Normal vectors of planar surface primitives do not depend on the point,
# but the normal method is generally called with a point.
normal(p::AbstractPlanarSurfacePrimitive, ::AbstractCoordinatePoint) = normal(p)
