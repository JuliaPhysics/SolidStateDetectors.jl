"""
    Box{T, CO} <: AbstractVolumePrimitive{T}

T: Type of values, e.g. Float64
CO: ClosedPrimitive or OpenPrimitive <-> whether surface belongs to it or not
"""
struct Box{T, CO} <: AbstractVolumePrimitive{T, CO}
    hX::T
    hY::T
    hZ::T
    origin::CartesianPoint{T}
    rotation::SMatrix{3,3,T,9}
end

Box{T, CO}(hX::Real, hY::Real, hZ::Real) where {T, CO<:Union{ClosedPrimitive, OpenPrimitive}} = 
    Box{T, CO}(T(hX), T(hY), T(hZ), CartesianPoint{T}(0,0,0), one(SMatrix{3, 3, T, 9}))
Box{T, CO}(hX::Real, hY::Real, hZ::Real, origin::CartesianPoint{T}) where {T, CO} = 
    Box{T, CO}(T(hX), T(hY), T(hZ), origin, one(SMatrix{3, 3, T, 9}))
   
ClosedPrimitive(b::Box{T}) where {T} = Box{T, ClosedPrimitive}(b.hX, b.hY, b.hZ, b.origin, b.rotation)
OpenPrimitive(b::Box{T}) where {T} = Box{T, OpenPrimitive}(b.hX, b.hY, b.hZ, b.origin, b.rotation)

extremum(b::Box{T}) where {T} = norm(CartesianPoint{T}(b.hX, b.hY, b.hZ))
rotation(p::AbstractVolumePrimitive) = p.rotation
origin(p::AbstractVolumePrimitive) = p.origin

scale(b::Box{T, CO}, s::SVector{3, <:Any}) where {T, CO} = Box{T, CO}( b.hX * s[1], b.hY * s[2], b.hZ * s[3], b.origin, b.rotation )
(*)(s::SVector{3, <:Any}, b::Box) = scale(b, s)

translate(b::Box{T, CO}, v::CartesianVector) where {T, CO} = Box{T, CO}(b.hX, b.hY, b.hZ, b.origin + v, b.rotation)
(+)(b::Box, v::CartesianVector) = translate(b, v)

rotate(b::Box{T, CO}, r::AbstractMatrix{T}) where {T, CO} = Box{T, CO}(b.hX, b.hY, b.hZ, r * b.origin, r * b.rotation)
# rotate_aroud_own_center(b::Box{T, CO}, r::AbstractMatrix{T}) where {T, CO} = Box{T, CO}(b.hX, b.hY, b.hZ, b.origin, r * b.rotation)
(*)(r::AbstractMatrix, b::Box) = rotate(b, r)

_in(pt::CartesianPoint, b::Box{<:Any, ClosedPrimitive}) =
    abs(pt.x) <= b.hX && abs(pt.y) <= b.hY && abs(pt.z) <= b.hZ
_in(pt::CartesianPoint, b::Box{<:Any, :OpenPrimitive}) = 
    abs(pt.x) < b.hX && abs(pt.y) < b.hY && abs(pt.z) < b.hZ
