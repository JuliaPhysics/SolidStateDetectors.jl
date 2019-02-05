# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).


struct CartesianPoint{ T <: RealQuantity } <: StaticArrays.FieldVector{3, T}
    x::T
    y::T
    z::T
end

struct CylindricalPoint{ T <: RealQuantity } <: StaticArrays.FieldVector{3, T}
    r::T
    θ::T # in radian
    z::T
end


function CylindricalPoint( pt::CartesianPoint{T} )::CylindricalPoint{T} where {T <: RealQuantity}
    cyl::CylindricalPoint{T} = CylindricalPoint(sqrt(pt.x * pt.x + pt.y * pt.y), atan(pt.y, pt.x), pt.z)
    while cyl.θ < 0 cyl = CylindricalPoint{T}(cyl.r, cyl.θ + 2π, cyl.z) end 
    while cyl.θ > 2π cyl = CylindricalPoint{T}(cyl.r, cyl.θ - 2π, cyl.z) end 
    while cyl.θ == 2π cyl = CylindricalPoint{T}(cyl.r, 0, cyl.z) end 
    return geom_round(cyl)
    # return cyl
end

function CartesianPoint( pt::CylindricalPoint{T} )::CartesianPoint{T} where {T <: RealQuantity}
    sθ::T, cθ::T = sincos(pt.θ)
    return geom_round(CartesianPoint{T}(pt.r * cθ, pt.r * sθ, pt.z))
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


function geom_round(pt::CylindricalPoint{T})::CylindricalPoint{T} where {T <: AbstractFloat} 
    return CylindricalPoint{T}( geom_round(pt.r), geom_round(pt.θ), geom_round(pt.z)  )
end

function geom_round(pt::CartesianPoint{T})::CartesianPoint{T} where {T <: AbstractFloat} 
    return CartesianPoint{T}( geom_round(pt.x), geom_round(pt.y), geom_round(pt.z)  )
end


# # Deprecated
# function CartFromCyl(p::CylindricalPoint{T}) where {T <: AbstractFloat}
#     CartesianPoint( p )
# end

# function CylFromCart(p::SVector{3, T})::CylindricalPoint{T} where {T <:AbstractFloat}
#     return CylindricalPoint(CartesianPoint{T}(p))
# end