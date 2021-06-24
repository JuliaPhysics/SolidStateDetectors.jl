@with_kw struct Torus{T,CO,TR,TP,TT} <: AbstractVolumePrimitive{T,CO}
    r_torus::T = 1
    r_tube::TR = 1 
    φ::TP = nothing
    θ::TT = nothing

    origin::CartesianPoint{T} = zero(CartesianPoint{T})
    rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})
end

Torus{T,CO,TR,TP,TT}( t::Torus{T,CO,TR,TP,TT}; COT = CO,
            origin::CartesianPoint{T} = t.origin,
            rotation::SMatrix{3,3,T,9} = t.rotation) where {T,CO<:Union{ClosedPrimitive, OpenPrimitive},TR,TP,TT} =
    Torus{T,COT,TR,TP,TT}(t.r_torus, t.r_tube, t.φ, t.θ, origin, rotation)

const FullTorus{T,CO} = Torus{T,CO,T,Nothing,Nothing}

function Geometry(::Type{T}, ::Type{Torus}, dict::AbstractDict, input_units::NamedTuple, transformations::Transformations{T}) where {T}
    length_unit = input_units.length
    angle_unit = input_units.angle
    r_torus = _parse_value(T, dict["r_torus"], length_unit)
    r_tube = _parse_radial_interval(T, dict["r_tube"], length_unit)
    φ = parse_φ_of_primitive(T, dict, angle_unit)
    θ = parse_θ_of_primitive(T, dict, angle_unit)
    t = Torus{T,ClosedPrimitive,typeof(r_tube),typeof(φ),typeof(θ)}(
        r_torus = r_torus, r_tube = r_tube, φ =φ, θ = θ)
    transform(t, transformations)
end

function surfaces(t::FullTorus{T,ClosedPrimitive}) where {T}
    tm = FullTorusMantle{T,:inwards}(t.r_torus, t.r_tube, t.φ, t.θ, t.origin, t.rotation)
    (tm, )
end
function surfaces(t::FullTorus{T,OpenPrimitive}) where {T}
    tm = FullTorusMantle{T,:outwards}(t.r_torus, t.r_tube, t.φ, t.θ, t.origin, t.rotation)
    (tm, )
end

# #Constructors
# function Torus(;r_torus = 1, r_tubeMin = 0, r_tubeMax = 1, φMin = 0, φMax = 2π, θMin = 0, θMax = 2π, z = 0)
#     T = float(promote_type(typeof.((r_torus, r_tubeMin, r_tubeMax, φMin, φMax, θMin, θMax, z))...))
#     r_tube = r_tubeMin == 0 ? T(r_tubeMax) : T(r_tubeMin)..T(r_tubeMax)
#     φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
#     θ = mod(T(θMax) - T(θMin), T(2π)) == 0 ? nothing : T(θMin)..T(θMax)
#     Torus( T, T(r_torus), r_tube, φ, θ, T(z))
# end

# Torus(r_torus, r_tubeMin, r_tubeMax, φMin, φMax, θMin, θMax, z) = Torus(;r_torus = r_torus, r_tubeMin = r_tubeMin, r_tubeMax = r_tubeMax, φMin = φMin, φMax = φMax, θMin = θMin, θMax = θMax, z = z)

# function Torus(r_torus::R1, r_tube::R2, z::TZ) where {R1<:Real, R2<:Real, TZ<:Real}
#     T = float(promote_type(R1, R2, TZ))
#     Torus( T, T(r_torus), T(r_tube), nothing, nothing, T(z))
# end

# function RoundChamfer(r_torus::R1, r_tube::R2, z::TZ) where {R1<:Real, R2<:Real, TZ<:Real}
#     T = float(promote_type(R1, R2, TZ))
#     Torus( T, T(r_torus), T(r_tube), nothing, T(0)..T(π/2), T(z))
# end

# in(p::AbstractCoordinatePoint, t::Torus{<:Any, <:Any, Nothing, Nothing}) =
#     _in_torr_r_tube(p, t.r_torus, t.r_tube, t.z)

# in(p::AbstractCoordinatePoint, t::Torus{<:Any, <:Any, <:AbstractInterval, Nothing}) =
#     _in_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_φ(p, t.φ)

# in(p::AbstractCoordinatePoint, t::Torus{<:Any, <:Any, Nothing, <:AbstractInterval}) =
#     _in_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_torr_θ(p, t.r_torus, t.θ, t.z)

# in(p::AbstractCoordinatePoint, t::Torus{<:Any, <:Any, <:AbstractInterval, <:AbstractInterval}) =
#     _in_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_φ(p, t.φ) && _in_torr_θ(p, t.r_torus, t.θ, t.z)
    

# # read-in
function Geometry(T::DataType, ::Type{Torus}, dict::AbstractDict, input_units::NamedTuple)
    length_unit = input_units.length
    angle_unit = input_units.angle
    r_torus::T = _parse_value(T, dict["r_torus"], length_unit)
    r_tube = _parse_radial_interval(T, dict["r_tube"], length_unit)
    φ = parse_φ_of_primitive(T, dict, angle_unit)
    θ = parse_θ_of_primitive(T, dict, angle_unit)
    z::T = ("z" in keys(dict) ? _parse_value(T, dict["z"], length_unit) : T(0))
    Torus(T, r_torus, r_tube, φ, θ, z)
end

# function Dictionary(t::Torus{T}) where {T}
#     dict = OrderedDict{String,Any}()
#     dict["r_torus"] = t.r_torus
#     dict["r_tube"] = typeof(t.r_tube) == T ? t.r_tube : OrderedDict{String,Any}("from" => t.r_tube.left, "to" => t.r_tube.right)
#     if !isnothing(t.φ) dict["phi"] = OrderedDict{String,Any}("from" => t.φ.left, "to" => t.φ.right) end
#     if !isnothing(t.θ) dict["theta"] = OrderedDict{String,Any}("from" => t.θ.left, "to" => t.θ.right) end
#     if t.z != 0 dict["z"] = t.z end
#     OrderedDict{String,Any}("torus" => dict)
# end


# get_r_tube_limits(t::Torus{T}) where {T} = (_left_radial_interval(t.r_tube),_right_radial_interval(t.r_tube))

# get_φ_limits(t::Torus{T, <:Any, Nothing, <:Any}) where {T} = (T(0), T(2π), true)
# get_φ_limits(t::Torus{T, <:Any, <:AbstractInterval, <:Any}) where {T} = (t.φ.left, t.φ.right, false)

# get_θ_limits(t::Torus{T, <:Any, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
# get_θ_limits(t::Torus{T, <:Any, <:Any, <:AbstractInterval}) where {T} = (t.θ.left, t.θ.right, false)

# function _is_torus_collapsed(t::Torus{T}) where {T}
#     r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
#     isapprox(r_tubeMin, r_tubeMax, atol = geom_atol_zero(T))
# end

# function _get_decomposed_surfaces_torus(t::Torus{T, T}) where {T}
#     r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
#     AbstractSurfacePrimitive[TorusMantle(t, r_tube = r_tubeMax)]
# end

# function _get_decomposed_surfaces_torus(t::Torus{T, <:AbstractInterval{T}}) where {T}
#     surfaces = AbstractSurfacePrimitive[]
#     r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
#     if isapprox(r_tubeMin, r_tubeMax, atol = geom_atol_zero(T))
#         push!(surfaces, TorusMantle(t, r_tube = r_tubeMax))
#     else
#         push!(surfaces, TorusMantle(t, r_tube = r_tubeMin), TorusMantle(t, r_tube = r_tubeMax))
#     end
#     surfaces
# end

# get_decomposed_surfaces(t::Torus{T, <:Any, Nothing, Nothing}) where {T} = _get_decomposed_surfaces_torus(t)

# function get_decomposed_surfaces(t::Torus{T, <:Any, <:AbstractInterval, Nothing}) where {T}
#     φMin::T, φMax::T, _ = get_φ_limits(t)
#     surfaces = _get_decomposed_surfaces_torus(t)
#     if !_is_torus_collapsed(t)
#         push!(surfaces, ToroidalAnnulus(t, φ = φMin), ToroidalAnnulus(t, φ = φMax))
#     end
#     surfaces
# end

# function get_decomposed_surfaces(t::Torus{T, <:Any, Nothing, <:AbstractInterval}) where {T}
#     r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
#     θMin::T, θMax::T, _ = get_θ_limits(t)
#     surfaces = _get_decomposed_surfaces_torus(t)
#     if !_is_torus_collapsed(t)
#         if r_tubeMin == T(0) && minmax(mod.((θMin, θMax), T(2π))...) == (T(0),T(π))
#             rMin = t.r_torus - r_tubeMax
#             rMax = t.r_torus + r_tubeMax
#             r = rMin == T(0) ? T(rMax) : T(rMin)..T(rMax)
#             push!(surfaces, CylindricalAnnulus(T, r, t.φ, t.z))
#         else
#             for θ in [θMin, θMax]
#                 mod(θ, T(2π)) in [T(0),T(π)] ? push!(surfaces, CylindricalAnnulus(t, θ = θ)) : push!(surfaces, ConeMantle(t, θ = θ))
#             end
#         end
#     end
#     surfaces
# end

# function get_decomposed_surfaces(t::Torus{T, <:Any, <:AbstractInterval, <:AbstractInterval}) where {T}
#     r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
#     θMin::T, θMax::T, _ = get_θ_limits(t)
#     φMin::T, φMax::T, _ = get_φ_limits(t)
#     surfaces = _get_decomposed_surfaces_torus(t)
#     if !_is_torus_collapsed(t)
#         if r_tubeMin == T(0) && minmax(mod.((θMin, θMax), T(2π))...) == (T(0),T(π))
#             rMin = t.r_torus - r_tubeMax
#             rMax = t.r_torus + r_tubeMax
#             r = rMin == T(0) ? T(rMax) : T(rMin)..T(rMax)
#             push!(surfaces, CylindricalAnnulus(T, r, t.φ, t.z))
#         else
#             for θ in [θMin, θMax]
#                 mod(θ, T(2π)) in [T(0),T(π)] ? push!(surfaces, CylindricalAnnulus(t, θ = θ)) : push!(surfaces, ConeMantle(t, θ = θ))
#             end
#         end
#         push!(surfaces, ToroidalAnnulus(t, φ = φMin), ToroidalAnnulus(t, φ = φMax))
#     end
#     surfaces
# end

# function sample(t::Torus{T}, step::Real) where {T}
#     r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
#     θMin::T, θMax::T, _ = get_θ_limits(t)
#     φMin::T, φMax::T, _ = get_φ_limits(t)
#     samples = [
#         CylindricalPoint{T}(t.r_torus+r_tube*cos(θ),φ,r_tube*sin(θ))
#         for r_tube in r_tubeMin:step:r_tubeMax
#         for θ in (r_tube == 0 ? θMin : θMin:step/r_tube:θMax)
#         for φ in (r_tube == 0 ? φMin : φMin:step/r_tube:φMax)
#     ]
# end

# function sample(t::Torus{T}, Nsamps::NTuple{3,Int} = (2,5,3)) where {T}
#     r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
#     θMin::T, θMax::T, _ = get_θ_limits(t)
#     φMin::T, φMax::T, _ = get_φ_limits(t)
#     samples = [
#         CylindricalPoint{T}(t.r_torus+r_tube*cos(θ),φ,r_tube*sin(θ))
#         for r_tube in (Nsamps[1] ≤ 1 ? r_tubeMin : range(r_tubeMin, r_tubeMax, length = Nsamps[1]))
#         for θ in (Nsamps[3] ≤ 1 ? θMin : range(θMin, θMax, length = Nsamps[3]))
#         for φ in (Nsamps[2] ≤ 1 ? φMin : range(φMin, φMax, length = Nsamps[2]))
#     ]
# end
