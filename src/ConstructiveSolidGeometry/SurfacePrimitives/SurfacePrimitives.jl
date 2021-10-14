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

function _get_n_points_in_arc_φ(p::AbstractSurfacePrimitive, n_arc::Int64)::Int64
    φMin, φMax = get_φ_limits(p)
    f = (φMax - φMin)/(2π)
    Int(ceil(n_arc*f))
end

function _get_n_points_in_arc_θ(p::AbstractSurfacePrimitive, n_arc::Int64)::Int64
    θMin, θMax = get_θ_limits(p)
    f = (θMax - θMin)/(2π)
    Int(ceil(n_arc*f))
end
