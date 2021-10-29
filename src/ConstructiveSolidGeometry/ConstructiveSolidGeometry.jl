# """
# # module ConstructiveSolidGeometry
# """
module ConstructiveSolidGeometry

    # Base Packages
    using LinearAlgebra

    # # Other Packages
    using IntervalSets
    using JSON
    using Parameters
    using PolygonOps
    using Polynomials
    using RecipesBase
    using Rotations
    using StaticArrays
    using Statistics
    using Unitful
    using YAML

    using DataStructures: OrderedDict

    import Base: in, *, +, -, &, size, zero

    abstract type AbstractCoordinateSystem end
    abstract type Cartesian <: AbstractCoordinateSystem end
    abstract type Cylindrical <: AbstractCoordinateSystem end
    const CoordinateSystemType = Union{Type{Cartesian}, Type{Cylindrical}}

    import Base: print, show
    print(io::IO, ::Type{Cartesian}) = print(io, "Cartesian")
    print(io::IO, ::Type{Cylindrical}) = print(io, "Cylindrical")
    show(io::IO, CS::CoordinateSystemType) = print(io, CS)
    show(io::IO,::MIME"text/plain", CS::CoordinateSystemType) = show(io, CS)

    # Tuples with ticks to sample with differently spaced ticks
    const CartesianTicksTuple{T} = NamedTuple{(:x,:y,:z), NTuple{3,Vector{T}}}
    const CylindricalTicksTuple{T} = NamedTuple{(:r,:φ,:z), NTuple{3,Vector{T}}}

    abstract type AbstractGeometry{T <: AbstractFloat} end

    abstract type AbstractPrimitive{T} <: AbstractGeometry{T} end
    abstract type ClosedPrimitive end
    abstract type OpenPrimitive end

    abstract type AbstractVolumePrimitive{T, CO} <: AbstractPrimitive{T} end
    abstract type AbstractSurfacePrimitive{T} <: AbstractPrimitive{T} end
    abstract type AbstractLinePrimitive{T} <: AbstractPrimitive{T} end

    abstract type AbstractConstructiveGeometry{T} <: AbstractGeometry{T} end

    include("Units.jl")
    include("PointsAndVectors/PointsAndVectors.jl")

    rotation(p::AbstractPrimitive) = p.rotation
    origin(p::AbstractPrimitive) = p.origin
    _transform_into_global_coordinate_system(pt::CartesianPoint, p::AbstractPrimitive) = (rotation(p) * pt) + origin(p)
    _transform_into_global_coordinate_system(pts::AbstractVector{<:CartesianPoint}, p::AbstractPrimitive) =
        broadcast(pt -> _transform_into_global_coordinate_system(pt, p), pts)
    _transform_into_object_coordinate_system(pt::CartesianPoint, p::AbstractPrimitive) = inv(rotation(p)) * (pt - origin(p)) 
    in(pt::CartesianPoint{T}, p::AbstractPrimitive{T}; csgtol::T = csg_default_tol(T)) where {T} = 
        _in(_transform_into_object_coordinate_system(pt, p), p; csgtol = csgtol)
    in(pt::CylindricalPoint{T}, p::AbstractPrimitive{T}; csgtol::T = csg_default_tol(T)) where {T} = 
        in(CartesianPoint(pt), p; csgtol = csgtol)

    """
        extreme_points(es::AbstractPrimitive{T}) where {T}

    Generic fallback of `extreme_points` for any primitive. 
    Returns 6 `CartesianPoint`s in both directions of each cartesian axes (x, y and z)
    around the origin of the primitive with distance determined by `extremum(es)`,
    which returns the maximum distance of an primitive to its origin.

    ## Arguments
    * `es::AbstractPrimitive{T}`: Any primitive, e.g. a `Box`.
    """
    function extreme_points(es::AbstractPrimitive{T}) where {T}
        o = origin(es)
        r = extremum(es)
        vX = r * CartesianVector{T}(one(T), zero(T), zero(T))
        vY = r * CartesianVector{T}(zero(T), one(T), zero(T))
        vZ = r * CartesianVector{T}(zero(T), zero(T), one(T))
        SVector{6,CartesianPoint{T}}(
            o - vX, o + vX, o - vY, o + vY, o - vZ, o + vZ,
        )
    end
    
    get_scale(es::AbstractPrimitive) = extremum(es)

    include("Transformations.jl")
    include("GeometryRounding.jl")
    include("VolumePrimitives/VolumePrimitives.jl")
    include("LinePrimitives/LinePrimitives.jl")
    include("SurfacePrimitives/SurfacePrimitives.jl")
    include("Intervals.jl")
    include("CSG.jl")
    include("IO.jl")

    # Plotting
    include("plotting/plotting.jl")

end
