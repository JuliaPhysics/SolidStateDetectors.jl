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

    import Base: in, *, +, -, &, size

    abstract type AbstractCoordinateSystem end
    abstract type Cartesian <: AbstractCoordinateSystem end
    abstract type Cylindrical <: AbstractCoordinateSystem end

    abstract type AbstractGeometry{T <: AbstractFloat} end

    abstract type AbstractPrimitive{T} <: AbstractGeometry{T} end
    abstract type AbstractVolumePrimitive{T} <: AbstractPrimitive{T} end
    abstract type AbstractSurfacePrimitive{T} <: AbstractPrimitive{T} end
    abstract type AbstractLinePrimitive{T} <: AbstractPrimitive{T} end

    abstract type AbstractConstructiveGeometry{T} <: AbstractGeometry{T} end

    include("Units.jl")
    include("PointsAndVectors.jl")
    include("GeometryRounding.jl")
    include("VolumePrimitives/VolumePrimitives.jl")
    include("SurfacePrimitives/SurfacePrimitives.jl")
    include("Sampling.jl")
    include("Transformations.jl")
    include("Intervals.jl")
    include("CSG.jl")
    include("IO.jl")
    include("Decompose.jl")

    # Plotting
    include("plotting/plotting.jl")

end
