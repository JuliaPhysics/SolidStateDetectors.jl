@with_kw struct RegularPrism{T,CO,N,TR} <: AbstractVolumePrimitive{T, CO}
    r::TR = 1
    hZ::T = 1

    origin::CartesianPoint{T} = zero(CartesianPoint{T})
    rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})
end

RegularPrism{T, CO, N, TR}( rp::RegularPrism{T, CO, N, TR}; COT = CO,
            origin::CartesianPoint{T} = rp.origin,
            rotation::SMatrix{3,3,T,9} = rp.rotation) where {T, CO<:Union{ClosedPrimitive, OpenPrimitive}, N, TR} =
    RegularPrism{T, COT, N, TR}(rp.r, rp.hZ, origin, rotation)

const TriangularPrism{T,CO,TR} = RegularPrism{T,CO,3,TR}
const QuadranglePrism{T,CO,TR} = RegularPrism{T,CO,4,TR}
const PentagonalPrism{T,CO,TR} = RegularPrism{T,CO,5,TR}
const HexagonalPrism{T,CO,TR}  = RegularPrism{T,CO,6,TR}

function Geometry(::Type{T}, ::Type{P}, dict::AbstractDict, input_units::NamedTuple, transformations::Transformations{T}
            ) where {T, P <: Union{TriangularPrism, QuadranglePrism, PentagonalPrism, HexagonalPrism}}
    length_unit = input_units.length
    angle_unit = input_units.angle
    origin = get_origin(T, dict, length_unit)
    rotation = get_rotation(T, dict, angle_unit)

    r = parse_r_of_primitive(T, dict, length_unit)
    hZ = if haskey(dict, "h")
        _parse_value(T, dict["h"], length_unit) / 2
    elseif haskey(dict, "z")
        z = parse_height_of_primitive(T, dict, length_unit)
        if z isa Real
            @warn "Deprecation warning: Field `z` for `RegularPrism` is deprecated. 
                Use `h` instead to specify the height of the primitive."
            z
        else # z isa Tuple
            @warn "Deprecation warning: Field `z` for `CRegularPrismone` is deprecated. 
                Use `h` instead to specify the height of the primitive.
                There might be a conflict with the possible field `origin`:
                The `z` component of the origin of the primitive is overwritten by the `z`."
            origin = CartesianPoint{T}(origin[1], origin[2], mean(z))
            (z[2] - z[1])/2
       end
    end
    
    g = if r isa Tuple # lazy workaround for now
        P{T,ClosedPrimitive,T}(r = r[2], hZ = hZ, origin = origin, rotation = rotation) -
        P{T,ClosedPrimitive,T}(r = r[1], hZ = T(1.1)*hZ, origin = origin, rotation = rotation)
    else
        P{T,ClosedPrimitive,T}(r = r, hZ = hZ, origin = origin, rotation = rotation)
    end
    transform(g, transformations)
end

const PrismAliases = Dict{Int, String}(
    3 => "TriangularPrism", 
    4 => "QuadranglePrism", 
    5 => "PentagonalPrism", 
    6 => "HexagonalPrism"
)

function Dictionary(rp::RegularPrism{T, <:Any, N})::OrderedDict{String, Any} where {T, N}
    dict = OrderedDict{String, Any}()
    dict["r"] = rp.r # always a Real 
    dict["h"] = rp.hZ*2
    if rp.origin != zero(CartesianVector{T}) dict["origin"] = rp.origin end
    if rp.rotation != one(SMatrix{3,3,T,9}) dict["rotation"] = Dictionary(rp.rotation) end
    OrderedDict{String, Any}(PrismAliases[N] => dict)
end

function vertices(rp::RegularPrism{T,ClosedPrimitive,N,T}) where {T,N}
    xys = [rp.r .* sincos(T(2π)*(n-1)/N) for n in 1:N]
    pts = [CartesianPoint{T}(xy[2], xy[1], z) for z in (-rp.hZ, rp.hZ) for xy in xys]
    _transform_into_global_coordinate_system(pts, rp)
end
function vertices(rp::RegularPrism{T,OpenPrimitive,N,T}) where {T,N}
    xys = [rp.r .* sincos(T(2π)*(n-1)/N) for n in N:-1:1]
    pts = [CartesianPoint{T}(xy[2], xy[1], z) for z in (-rp.hZ, rp.hZ) for xy in xys]
    _transform_into_global_coordinate_system(pts, rp)
end

function surfaces(rp::RegularPrism{T,<:Any,N,T}) where {T,N}
    vs = (vertices(rp))
    p_bot = Polygon{N,T}(vs[1:N])
    p_top = Polygon{N,T}(reverse(vs[N+1:end]))
    quads = Vector{Quadrangle{T}}(undef, N)
    for i in 1:N-1
        quads[i] = Quadrangle{T}((vs[i], vs[N+i], vs[N+i+1], vs[i+1]))
    end
    quads[N] = Quadrangle{T}((vs[N], vs[2N], vs[N+1], vs[1]))
    p_bot, p_top, quads...
end


function sample(rp::RegularPrism{T,<:Any,N,T})::Vector{CartesianPoint{T}} where {T,N}
    [vertices(rp)...]
end

function _in(pt::CartesianPoint{T}, rp::RegularPrism{T,ClosedPrimitive,N,T}; csgtol::T = csg_default_tol(T)) where {T,N} 
    abs(pt.z) <= rp.hZ + csgtol && begin
        r, φ = hypot(pt.x, pt.y), atan(pt.y, pt.x)
        α = T(π/N)
        _r = r * cos(α - mod(φ, 2α)) / cos(α) 
        _r <= rp.r + csgtol
    end
end
function _in(pt::CartesianPoint{T}, rp::RegularPrism{T,OpenPrimitive,N,T}; csgtol::T = csg_default_tol(T)) where {T,N} 
    abs(pt.z) < rp.hZ - csgtol && begin
        r, φ = hypot(pt.x, pt.y), atan(pt.y, pt.x)
        α = T(π/N)
        _r = r * cos(α - mod(φ, 2α)) / cos(α) 
        _r < rp.r - csgtol
    end
end

# # Convenience functions
# const TriangularPrism{T,TR,TZ} = RegularPrism{3,T,TR,TZ}
# const SquarePrism{T,TR,TZ}  = RegularPrism{4,T,TR,TZ}
# const PentagonalPrism{T,TR,TZ} = RegularPrism{5,T,TR,TZ}
# const HexagonalPrism{T,TR,TZ}  = RegularPrism{6,T,TR,TZ}

# TriangularPrism(args...) = RegularPrism(3, args...)
# SquarePrism(args...)  = RegularPrism(4, args...)
# PentagonalPrism(args...) = RegularPrism(5, args...)
# HexagonalPrism(args...)  = RegularPrism(6, args...)

# print(io::IO, rp::TriangularPrism{T, TR, TZ}) where {T,TR,TZ} = print(io, "TriangularPrism{$(T), $(TR), $(TZ)}($(rp.r), $(rp.z))")
# print(io::IO, rp::SquarePrism{T, TR, TZ})  where {T,TR,TZ} = print(io, "SquarePrism{$(T), $(TR), $(TZ)}($(rp.r), $(rp.z))")
# print(io::IO, rp::PentagonalPrism{T, TR, TZ}) where {T,TR,TZ} = print(io, "PentagonalPrism{$(T), $(TR), $(TZ)}($(rp.r), $(rp.z))")
# print(io::IO, rp::HexagonalPrism{T, TR, TZ})  where {T,TR,TZ} = print(io, "HexagonalPrism{$(T), $(TR), $(TZ)}($(rp.r), $(rp.z))")


# # Constructors
# function RegularPrism(N::Integer; rInner = 0, rOuter = 1, zMin = -1, zMax = 1)
#     T = float(promote_type(typeof.((rInner, rOuter, zMin, zMax))...))
#     r = rInner == 0 ? T(rOuter) : T(rInner)..T(rOuter)
#     z = zMax == -zMin ? T(zMax) : T(zMin)..T(zMax)
#     RegularPrism(N, T, r, z)
# end
# RegularPrism(N::Integer, rInner, rOuter, zMin, zMax) = RegularPrism(N; rInner = rInner, rOuter = rOuter, zMin = zMin, zMax = zMax)

# function RegularPrism(N::Integer, rOuter::R, height::H) where {R<:Real, H<:Real}
#     T = float(promote_type(R,H))
#     RegularPrism( N, T, T(rOuter), T(height/2))
# end


# @inline in(p::CylindricalPoint, hp::RegularPrism{N, T, <:Real}) where {N,T} = begin
#     _in_z(p, hp.z) && p.r * cos(T(π/N) - mod(p.φ, T(2π/N))) / cos(T(π/N)) <= hp.r
# end

# @inline in(p::CylindricalPoint, hp::RegularPrism{N, T, <:AbstractInterval}) where {N,T} = begin
#     _in_z(p, hp.z) && p.r * cos(T(π/N) - mod(p.φ, T(2π/N))) / cos(T(π/N)) in hp.r
# end

# @inline in(p::CartesianPoint, hp::RegularPrism) = in(CylindricalPoint(p), hp)

# # Special case: CartesianPoint in HexagonalPrism: use analytical formulas
# @inline in(p::CartesianPoint, hp::HexagonalPrism{T, <:Real}) where {T} =
#      _in_z(p, hp.z) && _in_y(p, hp.r * sqrt(T(3))/2) && _in_x(p, hp.r - abs(p.y)/sqrt(T(3)))

# @inline in(p::CartesianPoint, hp::HexagonalPrism{T, <:AbstractInterval{T}}) where {T} =
#     _in_z(p, hp.z) && abs(p.y) <= hp.r.right * sqrt(T(3))/2 &&
#     (
#         abs(p.y) >= hp.r.left * sqrt(T(3))/2 && _in_x(p, hp.r.right * T(0.5)) ||
#         abs(p.x) in (hp.r.left - abs(p.y) /sqrt(T(3)))..(hp.r.right - abs(p.y)/sqrt(T(3)))
#     )

# # read-in


# function Dictionary(rp::RegularPrism{N,T}) where {N,T}
#     dict = OrderedDict{String,Any}()
#     dict["r"] = typeof(rp.r) == T ? rp.r : OrderedDict{String,Any}("from" => rp.r.left, "to" => rp.r.right)
#     typeof(rp.z) == T ? dict["h"] = rp.z * 2 : dict["z"] = OrderedDict{String,Any}("from" => rp.z.left, "to" => rp.z.right)
#     OrderedDict{String,Any}(PrismAliases[N] => dict)
# end



# get_r_limits(rp::RegularPrism) = _radial_endpoints(rp.r)
# get_z_limits(rp::RegularPrism) = _linear_endpoints(rp.z)

# function get_decomposed_surfaces(rp::RegularPrism{N,T}) where {N, T}
#     rMin::T, rMax::T = get_r_limits(rp)
#     zMin::T, zMax::T = get_z_limits(rp)
#     surfaces = AbstractSurfacePrimitive[]
#     tol = geom_atol_zero(T)
#     if !isapprox(zMin, zMax, atol = tol)
#         if !isapprox(rMin, rMax, atol = tol)
#             if rMin == 0
#                 return AbstractSurfacePrimitive[RegularPolygon(rp, z = zMin), RegularPolygon(rp, z = zMax), RegularPrismMantle(rp, r = rMax)]
#             else
#                 return AbstractSurfacePrimitive[RegularPolygon(rp, z = zMin), RegularPolygon(rp, z = zMax), RegularPrismMantle(rp, r = rMin), RegularPrismMantle(rp, r = rMax)]
#             end
#         else
#             return AbstractSurfacePrimitive[RegularPrismMantle(rp, r = rMin)]
#         end
#     else
#         return AbstractSurfacePrimitive[RegularPolygon(rp, z = zMin)]
#     end
# end

# function sample(rp::RegularPrism{N,T}, Nsamps::NTuple{3,Int} = (2,N+1,2))::Vector{CylindricalPoint{T}} where {N,T}
#     rMin::T, rMax::T = get_r_limits(rp)
#     zMin::T, zMax::T = get_z_limits(rp)
#     samples = [
#         CylindricalPoint{T}(r*cos(π/N)/cos(π/N - mod(φ, 2π/N)),φ,z)
#         for z in (Nsamps[3] ≤ 1 ? zMin : range(zMin, zMax, length = Nsamps[3]))
#         for φ in (Nsamps[2] ≤ 1 ? 0 : range(0, 2π, length = Nsamps[2]))
#         for r in (Nsamps[1] ≤ 1 ? rMax : range(rMin, rMax, length = Nsamps[1]))
#     ]
# end
