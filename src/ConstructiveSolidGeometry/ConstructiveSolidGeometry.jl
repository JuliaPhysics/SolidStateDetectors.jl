"""
# module ConstructiveSolidGeometry
"""
module ConstructiveSolidGeometry

    # Base Packages
    using LinearAlgebra

    # # Other Packages
    using IntervalSets
    using Rotations
    using StaticArrays

    import Base: in, *
    
    abstract type AbstractCoordinateSystem end
    abstract type Cartesian <: AbstractCoordinateSystem end
    abstract type Cylindrical <: AbstractCoordinateSystem end

    abstract type AbstractGeometry{T} end 
    
    abstract type AbstractPrimitive{T} <: AbstractGeometry{T} end
    abstract type AbstractVolumePrimitive{T} <: AbstractPrimitive{T} end
    abstract type AbstractSurfacePrimitive{T} <: AbstractPrimitive{T} end
    abstract type AbstractLinePrimitive{T} <: AbstractPrimitive{T} end

    abstract type AbstractConstructiveGeometry{T} <: AbstractGeometry{T} end   

    include("PointsAndVectors.jl")
    include("GeometryRounding.jl")
    include("VolumePrimitives/VolumePrimitives.jl")
    include("Transformations.jl")
    include("Intervals.jl")
    include("CSG.jl")

end 
