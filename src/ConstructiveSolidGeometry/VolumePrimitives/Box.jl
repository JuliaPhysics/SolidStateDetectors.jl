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
 
function Geometry(::Type{T}, ::Type{Box}, dict::AbstractDict, input_units::NamedTuple) where {T}
    length_unit = input_units.length
    x = parse_interval_of_primitive(T, "x", dict, length_unit)
    y = parse_interval_of_primitive(T, "y", dict, length_unit)
    z = parse_interval_of_primitive(T, "z", dict, length_unit)
    μx = typeof(x) <: Real ? zero(T) : mean(x)
    μy = typeof(y) <: Real ? zero(T) : mean(y)
    μz = typeof(z) <: Real ? zero(T) : mean(z)
    origin = CartesianPoint{T}(μx, μy, μz)
    if haskey(dict, "translate")
        origin += parse_translate_vector(T, dict["translate"], length_unit)
    end
    scale = SVector{3, T}(1, 1, 1)
    if haskey(dict, "scale")
        # scale = parse_translate_vector(T, dict["translate"], length_unit)
        @info "Scale not yet implemented"
    end
    rot = if haskey(dict, "rotation")
        parse_rotation_matrix(T, dict["rotate"], length_unit)
    else
        one(SMatrix{3, 3, T, 9})
    end
    hX = typeof(x) <: Real ? x : width(x)/2
    hY = typeof(y) <: Real ? y : width(y)/2
    hZ = typeof(z) <: Real ? z : width(z)/2
    return Box{T, ClosedPrimitive}(
        scale[1] * hX, 
        scale[2] * hY, 
        scale[3] * hZ, 
        origin, rot
    )
end

function vertices(b::Box{T}) where {T}
    return SVector{8, CartesianPoint{T}}(
        b.rotation * SVector{3,T}(-b.hX, -b.hY, -b.hZ) .+ b.origin,
        b.rotation * SVector{3,T}(+b.hX, -b.hY, -b.hZ) .+ b.origin,
        b.rotation * SVector{3,T}(+b.hX, +b.hY, -b.hZ) .+ b.origin,
        b.rotation * SVector{3,T}(-b.hX, +b.hY, -b.hZ) .+ b.origin,
        b.rotation * SVector{3,T}(-b.hX, -b.hY, +b.hZ) .+ b.origin,
        b.rotation * SVector{3,T}(+b.hX, -b.hY, +b.hZ) .+ b.origin,
        b.rotation * SVector{3,T}(+b.hX, +b.hY, +b.hZ) .+ b.origin,
        b.rotation * SVector{3,T}(-b.hX, +b.hY, +b.hZ) .+ b.origin,
    )
end

sample(b::Box) = vertices(b)
