# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

abstract type AbstractCoordinatePoint{T, N, S} <: StaticArrays.FieldVector{N, T} end

struct CartesianPoint{ T <: RealQuantity } <: AbstractCoordinatePoint{T, 3, :Cartesian}
    x::T
    y::T
    z::T
end

struct CylindricalPoint{ T <: RealQuantity } <: AbstractCoordinatePoint{T, 3, :Cylindrical}
    r::T
    φ::T # in radian
    z::T
end


function CylindricalPoint( pt::CartesianPoint{T} )::CylindricalPoint{T} where {T <: RealQuantity}
    cyl::CylindricalPoint{T} = CylindricalPoint(sqrt(pt.x * pt.x + pt.y * pt.y), atan(pt.y, pt.x), pt.z)
    while cyl.φ < 0 cyl = CylindricalPoint{T}(cyl.r, cyl.φ + 2π, cyl.z) end # this should be removed at some point
    while cyl.φ > 2π cyl = CylindricalPoint{T}(cyl.r, cyl.φ - 2π, cyl.z) end # this should be removed at some point
    while cyl.φ == 2π cyl = CylindricalPoint{T}(cyl.r, 0, cyl.z) end # this should be removed at some point
    return cyl
end

function CartesianPoint( pt::CylindricalPoint{T} )::CartesianPoint{T} where {T <: RealQuantity}
    sφ::T, cφ::T = sincos(pt.φ)
    return CartesianPoint{T}(pt.r * cφ, pt.r * sφ, pt.z)
end



const SomeCartesianPoint{T} = Union{ # this should be removed
    CartesianPoint{T},
    StaticArray{Tuple{3},T},
}


const SomeCylindricalPoint = Union{ # this should be removed
    CylindricalPoint,
}


const AnyPoint{T} = Union{ # this should be removed -> AbstractCoordinatePoint{T}
    CartesianPoint{T},
    CylindricalPoint{T}
}
