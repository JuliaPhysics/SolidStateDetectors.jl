# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).


struct CartesianPoint{ T <: RealQuantity } <: StaticArrays.FieldVector{3, T}
    x::T
    y::T
    z::T
end


struct CylindricalPoint{ T <: RealQuantity } <: StaticArrays.FieldVector{3, T}
    r::T
    φ::T # in radian
    z::T
end


function CylindricalPoint( pt::CartesianPoint{T} )::CylindricalPoint{T} where {T <: RealQuantity}
    cyl::CylindricalPoint{T} = CylindricalPoint(sqrt(pt.x * pt.x + pt.y * pt.y), atan(pt.y, pt.x), pt.z)
    while cyl.φ < 0 cyl = CylindricalPoint{T}(cyl.r, cyl.φ + 2π, cyl.z) end
    while cyl.φ > 2π cyl = CylindricalPoint{T}(cyl.r, cyl.φ - 2π, cyl.z) end
    while cyl.φ == 2π cyl = CylindricalPoint{T}(cyl.r, 0, cyl.z) end
    return geom_round(cyl)
    # return cyl
end

function CartesianPoint( pt::CylindricalPoint{T} )::CartesianPoint{T} where {T <: RealQuantity}
    sφ::T, cφ::T = sincos(pt.φ)
    return geom_round(CartesianPoint{T}(pt.r * cφ, pt.r * sφ, pt.z))
end


const SomeCartesianPoint{T} = Union{
    CartesianPoint{T},
    StaticArray{Tuple{3},T},
}


const SomeCylindricalPoint = Union{
    CylindricalPoint,
}


const AnyPoint{T} = Union{
    CartesianPoint{T},
    CylindricalPoint{T}
}


# function geom_round(pt::CylindricalPoint{T})::CylindricalPoint{T} where {T <: AbstractFloat} 
#     return CylindricalPoint{T}( geom_round(pt.r), geom_round(pt.φ), geom_round(pt.z)  )
# end
#
# function geom_round(pt::CartesianPoint{T})::CartesianPoint{T} where {T <: AbstractFloat}
#     return CartesianPoint{T}( geom_round(pt.x), geom_round(pt.y), geom_round(pt.z)  )
# end


# # Deprecated
# function CartFromCyl(p::CylindricalPoint{T}) where {T <: AbstractFloat}
#     CartesianPoint( p )
# end

# function CylFromCart(p::SVector{3, T})::CylindricalPoint{T} where {T <:AbstractFloat}
#     return CylindricalPoint(CartesianPoint{T}(p))
# end
