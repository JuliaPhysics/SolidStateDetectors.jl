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

function _get_vert_lines_range(p::AbstractSurfacePrimitive{T}, n_arc::Int64, n_vert_lines::Int64)::Vector{Int64} where {T}
    φMin, φMax = get_φ_limits(p) 
    if mod(φMax-φMin, T(2π)) == 0 
        if n_vert_lines == 0 
            [] 
        elseif n_vert_lines == 1 
            [1] 
        else
            Int.(ceil.(range(1, (n_arc+1)*(n_vert_lines-1)/n_vert_lines, length = n_vert_lines)))
        end
    else
        Int.(ceil.(range(1, n_arc+1, length = max(n_vert_lines, 2))))
    end
end