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



# #Constructors
# function Sphere(; rMin = 0, rMax = 1)
#     T = float(promote_type(typeof.((rMin, rMax))...))
#     r = rMin == 0 ? T(rMax) : T(rMin)..T(rMax)
#     Sphere(T, r)
# end

# Sphere(rMin, rMax) = Sphere(; rMin = rMin, rMax = rMax)
# function Sphere(r::R) where {R <: Real}
#     T = float(R)
#     Sphere(T, T(r))
# end

# in(p::AbstractCoordinatePoint, s::Sphere) = _in_sph_r(p, s.r)

# # read-in
# function Geometry(::Type{T}, ::Type{Sphere}, dict::AbstractDict, input_units::NamedTuple) where {T}
#     length_unit = input_units.length
#     r = parse_r_of_primitive(T, dict, length_unit)
#     return Sphere(T, r)
# end

# function Dictionary(s::Sphere{T}) where {T}
#     dict = OrderedDict{String,Any}()
#     dict["r"] = typeof(s.r) == T ? s.r : OrderedDict{String,Any}("from" => s.r.left, "to" => s.r.right)
#     OrderedDict{String,Any}("sphere" => dict)
# end

# get_r_limits(s::Sphere{T}) where {T} = (_left_radial_interval(s.r), _right_radial_interval(s.r))

# get_decomposed_surfaces(s::Sphere{T, T}) where {T} = AbstractSurfacePrimitive[SphereMantle{T}(s.r)]

# function get_decomposed_surfaces(s::Sphere{T, <:AbstractInterval{T}}) where {T}
#     rMin::T, rMax::T = get_r_limits(s)
#     AbstractSurfacePrimitive[SphereMantle{T}(rMin), SphereMantle{T}(rMax)]
# end

# function sample(s::Sphere{T}, Nsamps::NTuple{3,Int} = (2,5,3))::Vector{CylindricalPoint{T}} where {T}
#     rMin::T, rMax::T = get_r_limits(s)
#     samples = [
#         CylindricalPoint{T}(r,φ,z)
#         for z in (Nsamps[3] ≤ 1 ? rMax : range(-rMax, rMax, length = Nsamps[3]))
#         for r in (Nsamps[1] ≤ 1 ? rMax : range((abs(z) > rMin ? 0 : sqrt(rMin^2 - z^2)), sqrt(rMax^2 - z^2), length = Nsamps[1]))
#         for φ in (Nsamps[2] ≤ 1 ? 0 : range(0, 2π, length = Nsamps[2]))
#     ]
# end
