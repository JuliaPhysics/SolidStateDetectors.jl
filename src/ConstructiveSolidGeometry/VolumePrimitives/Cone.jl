"""
    struct Cone{T,CO,RT,TP} <: AbstractVolumePrimitive{T, CO}

T: Type of values, e.g. Float64
CO: ClosedPrimitive or OpenPrimitive <-> whether surface belongs to it or not

* `r::TR`: 
    * TR = Real -> Solid Tube 
    * TR = (Real, Real) -> Pipe (r_in = r[1], r_out = r[2])
    * TR = ((Real,), (Real,)) Solid widening Tube  -> (r_bot = r[1][1], r_top = r[1][2])
    * TR = ((Real,Real), (Real,Real)) Solid widening Pipe ->\n(r_bot_in = r[1][1], r_bot_out = r[1][2], r_top_in = r[2][1], r_top_out = r[2][2])
    * TR = (Nothing, (Real,Real)) Cone ->\n(r_bot_in = r_bot_out = 0, r_top_in = r[2][1], r_top_out = r[2][2])
    * TR = ((Real,Real), Nothing) Cone ->\n(r_bot_in = r[1][1], r_bot_out = r[1][2], r_top_in = r_top_out = 0)
    * ... (+ elliptical cases -> (a, b))
    * Not all are implemented yet

* `φ::TP`: 
    * TP = Nothing <-> Full in φ
    * ...
* `zH::T`: half hight/length of the cone
"""
@with_kw struct Cone{T,CO,RT,TP} <: AbstractVolumePrimitive{T, CO}
    r::RT = 1
    φ::TP = nothing
    hZ::T = 1

    origin::CartesianPoint{T} = zero(CartesianPoint{T})
    rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})
end

Cone{T,CO,RT,TP}( c::Cone{T,CO,RT,TP}; COT = CO,
            origin::CartesianPoint{T} = b.origin,
            rotation::SMatrix{3,3,T,9} = b.rotation) where {T,CO<:Union{ClosedPrimitive, OpenPrimitive},RT,TP} =
    Cone{T,CO,RT,TP}(c.r, c.φ, c.hZ, origin, rotation)

const Tube{T,CO} = Cone{T,CO,T,Nothing}

_in(pt::CartesianPoint, c::Tube{<:Real, ClosedPrimitive}) = abs(pt.z) <= c.hZ && hypot(pt.x, pt.y) <= c.r
_in(pt::CartesianPoint, c::Tube{<:Real, OpenPrimitive})   = abs(pt.z) <  c.hZ && hypot(pt.x, pt.y) <  c.r


# #Constructors for Tubes
# Tube(;rMin = 0, rMax = 1, φMin = 0, φMax = 2π, zMin = -1/2, zMax = 1/2) = Cone(rMin, rMax, rMin, rMax, φMin, φMax, zMin, zMax)
# Tube(rMin, rMax, φMin, φMax, zMin, zMax) = Tube(; rMin = rMin, rMax = rMax, φMin = φMin, φMax = φMax, zMin = zMin, zMax = zMax)

# function Tube(r::R, height::H) where {R<:Real, H<:Real}
#     T = float(promote_type(R,H))
#     Cone(T, T(r), nothing, T(height)/2)
# end

# function Tube(rMin::R1, rMax::R2, height::H) where {R1<:Real, R2<:Real, H<:Real}
#     T = float(promote_type(R1,R2,H))
#     Cone(T, rMin == 0 ? T(rMax) : T(rMin)..T(rMax), nothing, T(height)/2)
# end

# # for Tubes
# get_r_at_z(c::Cone{T, <:Union{T, AbstractInterval{T}}, <:Any, <:Any}, z::Real) where {T} = c.r

# # for Cones
# get_r_at_z(c::Cone{T, Tuple{T,T}, <:Any, <:Any}, z::Real) where {T} = _get_r_at_z(c.r[1], c.r[2], c.z, z)

# function get_r_at_z(c::Cone{T, Tuple{I,I}, <:Any, <:Any}, z::Real) where {T, I<:AbstractInterval{T}}
#     r1::T = _get_r_at_z(c.r[1].left, c.r[2].left, c.z, z)
#     r2::T = _get_r_at_z(c.r[1].right, c.r[2].right, c.z, z)
#     r1..r2
# end

# function _get_r_at_z(rbot::TR, rtop::TR, cz::TZ, z::Real)::TR where {TR<:Real, TZ}
#     (rtop - rbot) * (z - _left_linear_interval(cz)) / _width_linear_interval(cz) + rbot
# end


# in(p::AbstractCoordinatePoint, c::Cone{<:Any, <:Any, Nothing, <:Any}) =
#     _in_z(p, c.z) && _in_cyl_r(p, get_r_at_z(c, p.z))

# in(p::AbstractCoordinatePoint, c::Cone{<:Any, <:Any, <:AbstractInterval, <:Any}) =
#     _in_z(p, c.z) && _in_φ(p, c.φ) && _in_cyl_r(p, get_r_at_z(c, p.z))

# # read-in
function Geometry(::Type{T}, t::Type{Cone}, dict::AbstractDict, input_units::NamedTuple, transformations::Transformations{T}) where {T}
    length_unit = input_units.length
    angle_unit = input_units.angle
    r = parse_r_of_primitive(T, dict, length_unit)
    φ = parse_φ_of_primitive(T, dict, angle_unit)
    h = parse_height_of_primitive(T, dict, length_unit)
    return Cone{T, ClosedPrimitive, typeof(r), typeof(φ)}(r = r, φ = φ, hZ = h/2)
end

function surfaces(t::Tube{T}) where {T}
    mantle = ConeMantle(t)
    e_bot = Ellipse{T,T,Nothing}(r = t.r, φ = nothing, origin = mantle.axis.origin, rotation = t.rotation)
    e_top = Ellipse{T,T,Nothing}(r = t.r, φ = nothing, origin = mantle.axis.origin + mantle.axis.direction, rotation = t.rotation)
    e_top, e_bot, mantle
end

# function Dictionary(g::Cone{T,<:Union{T, AbstractInterval}}) where {T}
#     dict = OrderedDict{String,Any}()
#     dict["r"] = typeof(g.r) == T ? g.r : OrderedDict{String,Any}("from" => g.r.left, "to" => g.r.right)
#     if !isnothing(g.φ) dict["phi"] = OrderedDict{String,Any}("from" => g.φ.left, "to" => g.φ.right) end
#     typeof(g.z) == T ? dict["h"] = g.z * 2 : dict["z"] = OrderedDict{String,Any}("from" => g.z.left, "to" => g.z.right) 
#     OrderedDict{String,Any}("tube" => dict)
# end

# function Dictionary(g::Cone{T,<:Tuple}) where {T}
#     dict = OrderedDict{String,Any}()
#     dict["r"] = OrderedDict{String,Any}()
#     dict["r"]["bottom"] = typeof(g.r[1]) == T ? g.r[1] : OrderedDict{String,Any}("from" => g.r[1].left, "to" => g.r[1].right)
#     dict["r"]["top"] = typeof(g.r[2]) == T ? g.r[2] : OrderedDict{String,Any}("from" => g.r[2].left, "to" => g.r[2].right)
#     if !isnothing(g.φ) dict["phi"] = OrderedDict{String,Any}("from" => g.φ.left, "to" => g.φ.right) end
#     typeof(g.z) == T ? dict["h"] = g.z * 2 : dict["z"] = OrderedDict{String,Any}("from" => g.z.left, "to" => g.z.right) 
#     OrderedDict{String,Any}("cone" => dict)
# end


# get_r_limits(c::Cone{T, <:Union{T, AbstractInterval{T}}, <:Any, <:Any}) where {T} =
#     (_left_radial_interval(c.r),_right_radial_interval(c.r),_left_radial_interval(c.r),_right_radial_interval(c.r))
# get_r_limits(c::Cone{T, <:Tuple, <:Any, <:Any}) where {T} =
#     (_left_radial_interval(c.r[1]),_right_radial_interval(c.r[1]),_left_radial_interval(c.r[2]),_right_radial_interval(c.r[2]))

# get_φ_limits(c::Cone{T, <:Any, Nothing, <:Any}) where {T} = (T(0), T(2π), true)
# get_φ_limits(c::Cone{T, <:Any, <:AbstractInterval, <:Any}) where {T} = (c.φ.left, c.φ.right, false)

# get_z_limits(c::Cone{T}) where {T} = (_left_linear_interval(c.z), _right_linear_interval(c.z))

# function _is_cone_collapsed(rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T, zMin::T, zMax::T) where {T}
#     tol = geom_atol_zero(T)
#     (isapprox(rbotMin, rbotMax, atol = tol) && isapprox(rtopMin, rtopMax, atol = tol)) || isapprox(zMin, zMax, atol = tol)
# end

# function _get_decomposed_surfaces_cone(c::Cone{T}, rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T, zMin::T, zMax::T) where {T}
#     surfaces = AbstractSurfacePrimitive[]
#     #top and bottom annulus
#     tol = geom_atol_zero(T)
#     if !isapprox(rbotMin, rbotMax, atol = tol)
#         push!(surfaces, CylindricalAnnulus(c, z = zMin))
#     end
#     if !isapprox(zMin, zMax, atol = tol)
#         if !isapprox(rtopMin, rtopMax, atol = tol)
#             push!(surfaces, CylindricalAnnulus(c, z = zMax))
#         end
#         #outer conemantle
#         push!(surfaces, ConeMantle(c, rbot = rbotMax, rtop = rtopMax))
#     end
#     surfaces
# end

# #2π Cones
# function get_decomposed_surfaces(c::Cone{T, <:Union{T, Tuple{T,T}}, Nothing, <:Any}) where {T}
#     rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = get_r_limits(c)
#     zMin::T, zMax::T = get_z_limits(c)
#     _get_decomposed_surfaces_cone(c, rbotMin, rbotMax, rtopMin, rtopMax, zMin, zMax)
# end

# function get_decomposed_surfaces(c::Cone{T, <:Union{<:AbstractInterval{T}, Tuple{I,I}}, Nothing, <:Any}) where {T, I<:AbstractInterval{T}}
#     rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = get_r_limits(c)
#     zMin::T, zMax::T = get_z_limits(c)
#     surfaces = _get_decomposed_surfaces_cone(c, rbotMin, rbotMax, rtopMin, rtopMax, zMin, zMax)
#     if !_is_cone_collapsed(rbotMin, rbotMax, rtopMin, rtopMax, zMin, zMax)
#         push!(surfaces, ConeMantle(c, rbot = rbotMin, rtop = rtopMin))
#     end
#     surfaces
# end

# #non 2π Cones
# function get_decomposed_surfaces(c::Cone{T, <:Union{T, Tuple{T,T}}, <:AbstractInterval{T}, <:Any}) where {T}
#     rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = get_r_limits(c)
#     zMin::T, zMax::T = get_z_limits(c)
#     φMin::T, φMax::T, _ = get_φ_limits(c)
#     surfaces = _get_decomposed_surfaces_cone(c, rbotMin, rbotMax, rtopMin, rtopMax, zMin, zMax)
#     if !_is_cone_collapsed(rbotMin, rbotMax, rtopMin, rtopMax, zMin, zMax)
#         push!(surfaces, ConalPlane(c, φ = φMin), ConalPlane(c, φ = φMax))
#     end
#     surfaces
# end

# function get_decomposed_surfaces(c::Cone{T, <:Union{<:AbstractInterval{T}, Tuple{I,I}}, <:AbstractInterval{T}, <:Any}) where {T, I<:AbstractInterval{T}}
#     rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = get_r_limits(c)
#     zMin::T, zMax::T = get_z_limits(c)
#     φMin::T, φMax::T, _ = get_φ_limits(c)
#     surfaces = _get_decomposed_surfaces_cone(c, rbotMin, rbotMax, rtopMin, rtopMax, zMin, zMax)
#     if !_is_cone_collapsed(rbotMin, rbotMax, rtopMin, rtopMax, zMin, zMax)
#         push!(surfaces, ConalPlane(c, φ = φMin), ConalPlane(c, φ = φMax))
#         push!(surfaces, ConeMantle(c, rbot = rbotMin, rtop = rtopMin))
#     end
#     surfaces
# end

# function sample(c::Cone{T}, step::Real)::Vector{CylindricalPoint{T}} where {T}
#     zMin::T, zMax::T = get_z_limits(c)
#     φMin::T, φMax::T, _ = get_φ_limits(c)
#     samples = [
#         CylindricalPoint{T}(r,φ,z)
#         for z in zMin:step:zMax
#         for r in _left_radial_interval(get_r_at_z(c, z)):step:_right_radial_interval(get_r_at_z(c, z))
#         for φ in (r == 0 ? φMin : φMin:step/r:φMax)
#     ]
# end

# function sample(c::Cone{T}, Nsamps::NTuple{3,Int})::Vector{CylindricalPoint{T}} where {T}
#     zMin::T, zMax::T = get_z_limits(c)
#     φMin::T, φMax::T, _ = get_φ_limits(c)
#     samples = [
#         CylindricalPoint{T}(r,φ,z)
#         for z in (Nsamps[3] ≤ 1 ? zMin : range(zMin, zMax, length = Nsamps[3]))
#         for r in (Nsamps[1] ≤ 1 ? _left_radial_interval(get_r_at_z(c, z)) : range(_left_radial_interval(get_r_at_z(c, z)), _right_radial_interval(get_r_at_z(c, z)), length = Nsamps[1]))
#         for φ in (Nsamps[2] ≤ 1 ? φMin : range(φMin, φMax, length = Nsamps[2]))
#     ]
# end

# @inline sample(c::Cone{T}) where {T} = sample(c, (2,3,3))
# @inline sample(c::Cone{T, <:Any, Nothing}) where {T} = sample(c, (2,5,3))
