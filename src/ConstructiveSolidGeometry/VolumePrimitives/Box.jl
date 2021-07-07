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

Box{T, CO}( b::Box{T, CO}; COT = CO,
            origin::CartesianPoint{T} = b.origin,
            rotation::SMatrix{3,3,T,9} = b.rotation) where {T, CO<:Union{ClosedPrimitive, OpenPrimitive}} =
    Box{T, COT}(b.hX, b.hY, b.hZ, origin, rotation)

function _in(pt::CartesianPoint{T}, b::Box{<:Any, ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T}
    abs(pt.x) <= b.hX + csgtol && 
    abs(pt.y) <= b.hY + csgtol && 
    abs(pt.z) <= b.hZ + csgtol 
end

_in(pt::CartesianPoint{T}, b::Box{<:Any, OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} = 
    abs(pt.x) < b.hX - csgtol && abs(pt.y) < b.hY - csgtol && abs(pt.z) < b.hZ - csgtol
 

function Geometry(::Type{T}, ::Type{Box}, dict::AbstractDict, input_units::NamedTuple, transformations::Transformations{T}) where {T}
    length_unit = input_units.length
    origin = get_origin(T, dict, length_unit)
    rotation = get_rotation(T, dict, input_units.angle) 
    hX, hY, hZ = if haskey(dict, "widths")
        @assert length(dict["widths"]) == 3
        _parse_value(T, dict["widths"], length_unit)/2
    elseif haskey(dict, "halfwidths")
        @assert length(dict["halfwidths"]) == 3
        _parse_value(T, dict["halfwidths"], length_unit)
    elseif haskey(dict, "hX") && haskey(dict, "hZ") && haskey(dict, "hZ")
        _parse_value(T, dict["hX"], length_unit), 
        _parse_value(T, dict["hY"], length_unit), 
        _parse_value(T, dict["hZ"], length_unit)
    elseif haskey(dict, "x") && haskey(dict, "y") && haskey(dict, "z")
        @warn "Deprecation warning: Detected old primitive definition for `Box`. 
            Please update your configuration file to the new format 
            via `widths`, `halfwidths` or `hX`, `hY` and `hZ`
            in combination with possible fields `origin` and `rotate`.
            The old definition overwrites the optional field `origin`."        
        x = parse_interval_of_primitive(T, "x", dict, length_unit)
        y = parse_interval_of_primitive(T, "y", dict, length_unit)
        z = parse_interval_of_primitive(T, "z", dict, length_unit)
        μx = typeof(x) <: Real ? zero(T) : mean(x)
        μy = typeof(y) <: Real ? zero(T) : mean(y)
        μz = typeof(z) <: Real ? zero(T) : mean(z)
        origin = CartesianPoint{T}(μx, μy, μz)
        hX = typeof(x) <: Real ? x : (x[2] - x[1])/2
        hY = typeof(y) <: Real ? y : (y[2] - y[1])/2
        hZ = typeof(z) <: Real ? z : (z[2] - z[1])/2
        hX, hY, hZ
    end
    box = Box{T, ClosedPrimitive}(
        hX = hX, 
        hY = hY, 
        hZ = hZ, 
        origin = origin,
        rotation = rotation
    )
    transform(box, transformations)
end

function vertices(b::Box{T}) where {T}
    return (
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

function sample(b::Box{T})::Vector{CartesianPoint{T}} where {T} 
    [vertices(b)...]
end

function surfaces(b::Box{T}) where {T}
    vs = vertices(b)
    return (
        Quadrangle{T}((vs[1], vs[2], vs[3], vs[4])),
        Quadrangle{T}((vs[5], vs[6], vs[2], vs[1])),
        Quadrangle{T}((vs[8], vs[7], vs[6], vs[5])),
        Quadrangle{T}((vs[6], vs[7], vs[3], vs[2])),
        Quadrangle{T}((vs[7], vs[8], vs[4], vs[3])),
        Quadrangle{T}((vs[8], vs[5], vs[1], vs[4])),
    )
end

