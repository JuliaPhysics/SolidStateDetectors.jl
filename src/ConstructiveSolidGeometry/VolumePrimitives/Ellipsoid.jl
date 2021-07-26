@with_kw struct Ellipsoid{T,CO,TR,TP,TT} <: AbstractVolumePrimitive{T, CO}
    r::TR = 1
    φ::TP = nothing
    θ::TT = nothing

    origin::CartesianPoint{T} = zero(CartesianPoint{T})
    rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})
end

Ellipsoid{T,CO,TR,TP,TT}( e::Ellipsoid{T,CO,TR,TP,TT}; COT = CO,
            origin::CartesianPoint{T} = e.origin,
            rotation::SMatrix{3,3,T,9} = e.rotation) where {T,CO<:Union{ClosedPrimitive, OpenPrimitive},TR,TP,TT} =
    Ellipsoid{T,COT,TR,TP,TT}(e.r, e.φ, e.θ, origin, rotation)

const Sphere{T,CO,TP,TT} = Ellipsoid{T,CO,T,TP,TT}
const FullSphere{T,CO} = Ellipsoid{T,CO,T,Nothing,Nothing}
const FullEllipsoid{T,CO} = Ellipsoid{T,CO,NTuple{3,T},Nothing,Nothing}

function Geometry(::Type{T}, ::Type{Ellipsoid}, dict::AbstractDict, input_units::NamedTuple, transformations::Transformations{T}) where {T}
    length_unit = input_units.length
    angle_unit = input_units.angle
    origin = get_origin(T, dict, length_unit)
    rotation = get_rotation(T, dict, angle_unit)
    
    r = parse_r_of_primitive(T, dict, input_units.length)
    φ = parse_φ_of_primitive(T, dict, angle_unit)
    if !(φ === nothing)
        error("Partial Ellipsoid (`φ = φ`) is not yet supported.")
    end
    θ = parse_θ_of_primitive(T, dict, angle_unit)
    if !(θ === nothing)
        error("Partial Ellipsoid (`θ = θ`) is not yet supported.")
    end
    e = if r isa Tuple{T,T}
        Ellipsoid{T,ClosedPrimitive,T,typeof(φ),typeof(θ)}(
            r = r[2], 
            φ = φ, 
            θ = θ,
            origin = origin,
            rotation = rotation
        ) - Ellipsoid{T,OpenPrimitive,T,typeof(φ),typeof(θ)}(
            r = r[1], 
            φ = φ, 
            θ = θ,
            origin = origin,
            rotation = rotation
        )
    else
        Ellipsoid{T,ClosedPrimitive,T,typeof(φ),typeof(θ)}(
            r = r, 
            φ = φ, 
            θ = θ,
            origin = origin,
            rotation = rotation
        )
    end
    transform(e, transformations)
end

function Dictionary(e::Ellipsoid{T})::OrderedDict{String, Any} where {T}
    dict = OrderedDict{String, Any}()
    dict["r"] = e.r # always a Real 
    if !isnothing(e.φ) error("Partial Ellipsoid (`φ = φ`) is not yet supported.") end
    if !isnothing(e.θ) error("Partial Ellipsoid (`θ = θ`) is not yet supported.") end
    if e.origin != zero(CartesianVector{T}) dict["origin"] = e.origin end
    if e.rotation != one(SMatrix{3,3,T,9}) dict["rotation"] = Dictionary(e.rotation) end
    OrderedDict{String, Any}("sphere" => dict)
end

_in(pt::CartesianPoint{T}, s::FullSphere{<:Any, ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T} =
    hypot(pt.x, pt.y, pt.z) <= s.r + csgtol

_in(pt::CartesianPoint{T}, s::FullSphere{<:Any, OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T} =
    hypot(pt.x, pt.y, pt.z) < s.r - csgtol

function surfaces(e::Ellipsoid{T,ClosedPrimitive}) where {T}
    em = EllipsoidMantle{T,typeof(e.r),typeof(e.φ),typeof(e.θ),:inwards}(e.r, e.φ, e.θ, e.origin, e.rotation)
    (em,)
end
function surfaces(e::Ellipsoid{T,OpenPrimitive}) where {T}
    em = EllipsoidMantle{T,typeof(e.r),typeof(e.φ),typeof(e.θ),:outwards}(e.r, e.φ, e.θ, e.origin, e.rotation)
    (em,)
end
