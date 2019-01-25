# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).


struct CartesianPoint{T<:RealQuantity} <: StaticArrays.FieldVector{3,T}
    x::T
    y::T
    z::T
end



struct CylindricalPoint{T<:Real}
    r::T
    θ::T # in radian
    z::T
end



const SomeCartesianPoint{T} = Union{
    CartesianPoint{T},
    StaticArray{Tuple{3},T},
}


const SomeCylindricalPoint = Union{
    CylindricalPoint,
    CoordinateTransformations.Cylindrical
}


const SpacialPoint = Union{SomeCartesianPoint, SomeCylindricalPoint}



function CylindricalPoint{T}(p::StaticArray{Tuple{3}}) where {T}
    x = p[1]; y = p[2]; z = p[3]
    r = sqrt(x * x + y * y)
    θ = atan(y, x)
    CylindricalPoint{T}(r, θ, z)
end


@inline CylindricalPoint(p::StaticArray{Tuple{3},T,1}) where {T<:Real} =
    CylindricalPoint{float(T)}(p)

@inline Base.convert(T::Type{<:CylindricalPoint}, p::StaticArray{Tuple{3}}) = T(p)


@inline CartesianPoint{T}(p::CylindricalPoint) where {T} = CartesianPoint{T}(p.r * cos(p.θ), p.r * sin(p.θ), p.z)
@inline SVector{3,T}(p::CylindricalPoint) where {T} = SVector(CartesianPoint{T}(p))

@inline CartesianPoint(p::CylindricalPoint{T}) where {T} = CartesianPoint{T}(p)
@inline SVector{3}(p::CylindricalPoint{T}) where {T} = SVector{3,float(T)}(p)

@inline Base.convert(T::Type{<:StaticArray{Tuple{3}}}, p::CylindricalPoint) = T(p)



@inline CylindricalPoint{T}(p::CoordinateTransformations.Cylindrical) where {T} =
    CylindricalPoint{T}(p.r, p.θ, p.z)

@inline CylindricalPoint(p::CoordinateTransformations.Cylindrical{T}) where {T} = CylindricalPoint{T}(p)

@inline Base.convert(T::Type{<:CylindricalPoint}, p::CoordinateTransformations.Cylindrical) = T(p)


CoordinateTransformations.Cylindrical{T}(p::CylindricalPoint) where {T} =
    CoordinateTransformations.Cylindrical(p.r, p.θ, p.z)

CoordinateTransformations.Cylindrical(p::CylindricalPoint{T}) where {T} =
    CoordinateTransformations.Cylindrical{T}(p)

@inline Base.convert(T::Type{<:CoordinateTransformations.Cylindrical}, p::CylindricalPoint) = T(p)
