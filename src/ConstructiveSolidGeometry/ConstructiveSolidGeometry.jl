"""
# module ConstructiveSolidGeometry
"""
module ConstructiveSolidGeometry

    # Base Packages
    using LinearAlgebra

    # # Other Packages
    using IntervalSets
    using RecipesBase
    using Rotations
    using StaticArrays
    using Unitful
    using Statistics

    using DataStructures: OrderedDict

    import Base: in, *, +, -, &, size

    abstract type AbstractCoordinateSystem end
    abstract type Cartesian <: AbstractCoordinateSystem end
    abstract type Cylindrical <: AbstractCoordinateSystem end
    const CoordinateSystemType = Union{Type{Cartesian}, Type{Cylindrical}}

    # import Base: print, show
    # print(io::IO, ::Type{Cartesian}) = print(io, "Cartesian")
    # print(io::IO, ::Type{Cylindrical}) = print(io, "Cylindrical")
    # show(io::IO, CS::CoordinateSystemType) = print(io, CS)
    # show(io::IO,::MIME"text/plain", CS::CoordinateSystemType) = show(io, CS)

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
    abstract type AbstractTransformedGeometry{T} <: AbstractGeometry{T} end

    include("Units.jl")
    include("PointsAndVectors/PointsAndVectors.jl")

    _transform_into_object_coordinate_system(pt::CartesianPoint, p::AbstractVolumePrimitive) = inv(rotation(p)) * (pt - origin(p)) 
    in(pt::CartesianPoint, p::AbstractVolumePrimitive) = _in(_transform_into_object_coordinate_system(pt, p), p)
    # Do we want to store the rotation matrix permanently in the primitive?
    # We should do tests regarding the performance. It can be easily added later. 

    include("GeometryRounding.jl")
    include("VolumePrimitives/VolumePrimitives.jl")
    include("LinePrimitives/LinePrimitives.jl")
    include("SurfacePrimitives/SurfacePrimitives.jl")
    include("Transformations.jl")
    include("Intervals.jl")
    include("CSG.jl")
    include("IO.jl")
    include("Sampling.jl")

    # Plotting
    include("plotting/plotting.jl")

end
