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
    r = parse_r_of_primitive(T, dict, length_unit)
    z = parse_height_of_primitive(T, dict, length_unit)
    hZ = typeof(z) <: Real ? z : width(z)/2
    origin = if typeof(z) <: Real
        CartesianPoint{T}(zero(T), zero(T), zero(T))
    else
        CartesianPoint{T}(zero(T), zero(T), mean(z))
    end
    return P{T,ClosedPrimitive,typeof(r)}(r = r, hZ = hZ, origin = origin)
end


# const PrismAliases = Dict{Int, String}(
#     3 => "TriangularPrism", 
#     4 => "SquarePrism", 
#     5 => "PentagonalPrism", 
#     6 => "HexagonalPrism"
# )

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



# get_r_limits(rp::RegularPrism) = (_left_radial_interval(rp.r), _right_radial_interval(rp.r))
# get_z_limits(rp::RegularPrism) = (_left_linear_interval(rp.z), _right_linear_interval(rp.z))

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
