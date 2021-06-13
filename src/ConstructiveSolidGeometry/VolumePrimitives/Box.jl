"""
    Box{T, CO} <: AbstractVolumePrimitive{T}

T: Type of values, e.g. Float64
CO: ClosedPrimitive or OpenPrimitive <-> whether surface belongs to it or not
"""
@with_kw struct Box{T, CO} <: AbstractVolumePrimitive{T, CO}
    hX::T = 1
    hY::T = 1
    hZ::T = 1
    origin::CartesianPoint{T} = zero(CartesianPoint{T})
    rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})
end

Box{T, CO}( b::Box{T, CO}; 
            COT = CO,
            scale::SVector{3,T} = ones(SVector{3,T}),
            origin::CartesianPoint{T} = zero(CartesianPoint{T}),
            rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})) where {T, CO<:Union{ClosedPrimitive, OpenPrimitive}} = 
    Box{T, COT}(b.hX * scale[1], b.hY * scale[2], b.hZ * scale[3], origin, rotation)

extremum(b::Box{T}) where {T} = norm(CartesianPoint{T}(b.hX, b.hY, b.hZ))

_in(pt::CartesianPoint, b::Box{<:Any, ClosedPrimitive}) =
    abs(pt.x) <= b.hX && abs(pt.y) <= b.hY && abs(pt.z) <= b.hZ
_in(pt::CartesianPoint, b::Box{<:Any, :OpenPrimitive}) = 
    abs(pt.x) < b.hX && abs(pt.y) < b.hY && abs(pt.z) < b.hZ
 
function Geometry(::Type{T}, ::Type{Box}, dict::AbstractDict, input_units::NamedTuple, transformations) where {T}
    length_unit = input_units.length
    x = parse_interval_of_primitive(T, "x", dict, length_unit)
    y = parse_interval_of_primitive(T, "y", dict, length_unit)
    z = parse_interval_of_primitive(T, "z", dict, length_unit)
    μx = typeof(x) <: Real ? zero(T) : mean(x)
    μy = typeof(y) <: Real ? zero(T) : mean(y)
    μz = typeof(z) <: Real ? zero(T) : mean(z)
    origin = CartesianPoint{T}(μx, μy, μz)
    scale = ones(SVector{3,T})
    hX = typeof(x) <: Real ? x : width(x)/2
    hY = typeof(y) <: Real ? y : width(y)/2
    hZ = typeof(z) <: Real ? z : width(z)/2
    box = Box{T, ClosedPrimitive}(
        hX = scale[1] * hX, 
        hY = scale[2] * hY, 
        hZ = scale[3] * hZ, 
        origin = origin
    )
    transform(box, transformations)
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

function surfaces(b::Box{T}) where {T}
    vs = vertices(b)
    return SVector{6, Quadrangle{T}}(
        Quadrangle{T}((vs[1], vs[2], vs[3], vs[4])),
        Quadrangle{T}((vs[5], vs[6], vs[2], vs[1])),
        Quadrangle{T}((vs[8], vs[7], vs[6], vs[5])),
        Quadrangle{T}((vs[6], vs[7], vs[3], vs[2])),
        Quadrangle{T}((vs[7], vs[8], vs[4], vs[3])),
        Quadrangle{T}((vs[8], vs[5], vs[1], vs[4])),
    )
end

