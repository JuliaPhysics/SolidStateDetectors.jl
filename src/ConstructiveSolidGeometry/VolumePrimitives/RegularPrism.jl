"""
    struct RegularPrism{T,CO,N,TR} <: AbstractVolumePrimitive{T, CO}

Volume primitive describing a [Prism](@ref) with base plates are regular polygons
which are parallel to the `xy` plane. If the regular polygon base plate is projected to 
the `xy` plane, one of the vertices lays on the `x` axis.


## Parametric types
* `T`: Precision type.
* `CO`: Describes whether the surface belongs to the primitive. 
    It can be `ClosedPrimitive`, i.e. the surface points belong to the primitive,
    or `OpenPrimitive`, i.e. the surface points do not belong to the primitive.
* `N`: Number of vertices of the regular polygon that defines the base of the prism.
* `TR`: Type of `r`.
    * `TR == T`: Regular polygon base (all vertices have the same distance to the center).
    
## Fields
* `r::TR`: Distance of the vertices to the center of the regular polygon base (in m).
* `hZ::T`: Half of the width in `z` dimension (in m).
* `origin::CartesianPoint{T}`: The position of the center of the `RegularPrism`.
* `rotation::SMatrix{3,3,T,9}`: Matrix that describes a rotation of the `RegularPrism` around its `origin`.

## Definition in Configuration File

So far, only `HexagonalPrism` can be defined in the configuration files.
A `HexagonalPrism` is defined in the configuration file as part of the `geometry` field 
of an object through the field `HexagonalPrism`.

An example definition of a `HexagonalPrism` looks like this:
```yaml
HexagonalPrism:
  r: 1.0 # => r = 1.0
  h: 2.0 # => hZ = 1.0
```

See also [Constructive Solid Geometry (CSG)](@ref).
"""
struct RegularPrism{N,T,CO,TR} <: AbstractVolumePrimitive{T, CO}
    r::TR
    hZ::T

    origin::CartesianPoint{T}
    rotation::SMatrix{3,3,T,9}
end

#Type conversion happens here
function RegularPrism{N,T,CO}(r, hZ, origin, rotation) where {T,CO,N}
    _r = _csg_convert_args(T, r)
    _hZ = _csg_convert_args(T, hZ)
    RegularPrism{N,T,CO,typeof(_r)}(_r, _hZ, origin, rotation)
end

#Type promotion happens here
function RegularPrism{N}(CO, r::TR, hZ::TZ, origin::PT, rotation::ROT) where {N,TR, TZ, PT, ROT}
    eltypes = _csg_get_promoted_eltype.((TR, TZ, PT, ROT))
    T = float(promote_type(eltypes...))
    RegularPrism{N,T,CO}(r, hZ, origin, rotation)
end

function RegularPrism{N}(::Type{CO}=ClosedPrimitive;
    r = 1, 
    hZ = 1,
    origin = zero(CartesianPoint{Int}), 
    rotation = one(SMatrix{3, 3, Int, 9})
) where {N, CO<:Union{ClosedPrimitive,OpenPrimitive}}
    RegularPrism{N}(CO, r, hZ, origin, rotation)
end

function RegularPrism{N,T}(::Type{CO}=ClosedPrimitive;
    r = 1.0, 
    hZ = 1.0,
    origin = zero(CartesianPoint{Float64}), 
    rotation = one(SMatrix{3, 3, Float64, 9})
) where {N,T <:Real, CO}
    RegularPrism{N,T,CO}(r, hZ, origin, rotation)
end

RegularPrism{N,T, CO, TR}( rp::RegularPrism{N,T, CO, TR}; COT = CO,
            origin::CartesianPoint{T} = rp.origin,
            rotation::SMatrix{3,3,T,9} = rp.rotation) where {T, CO<:Union{ClosedPrimitive, OpenPrimitive}, N, TR} =
    RegularPrism{N,T, COT, TR}(rp.r, rp.hZ, origin, rotation)

const TriangularPrism{T,CO,TR} = RegularPrism{3,T,CO,TR}
const QuadranglePrism{T,CO,TR} = RegularPrism{4,T,CO,TR}
const PentagonalPrism{T,CO,TR} = RegularPrism{5,T,CO,TR}
const HexagonalPrism{T,CO,TR}  = RegularPrism{6,T,CO,TR}

_get_N_prism(::Type{T},::Type{TriangularPrism},CO,r,hZ,origin,rotation) where {T} = RegularPrism{3,T}(CO, r = r, hZ = hZ, origin = origin, rotation = rotation)
_get_N_prism(::Type{T},::Type{QuadranglePrism},CO,r,hZ,origin,rotation) where {T} = RegularPrism{4,T}(CO, r = r, hZ = hZ, origin = origin, rotation = rotation)
_get_N_prism(::Type{T},::Type{PentagonalPrism},CO,r,hZ,origin,rotation) where {T} = RegularPrism{5,T}(CO, r = r, hZ = hZ, origin = origin, rotation = rotation)
_get_N_prism(::Type{T},::Type{HexagonalPrism},CO,r,hZ,origin,rotation) where {T} = RegularPrism{6,T}(CO, r = r, hZ = hZ, origin = origin, rotation = rotation)
            
function Geometry(::Type{T}, ::Type{P}, dict::AbstractDict, input_units::NamedTuple, transformations::Transformations{T}
            ) where {T, P <: Union{TriangularPrism, QuadranglePrism, PentagonalPrism, HexagonalPrism}}
    length_unit = input_units.length
    angle_unit = input_units.angle
    origin = get_origin(T, dict, length_unit)
    rotation = get_rotation(T, dict, angle_unit)

    r = parse_r_of_primitive(T, dict, length_unit)
    @assert haskey(dict,"h") || haskey(dict,"z") "Please specify 'h' or 'z'."
    hZ = if haskey(dict, "h")
        _parse_value(T, dict["h"], length_unit) / 2
    end
    
    g = if r isa Tuple # lazy workaround for now
        _get_N_prism(T,P,ClosedPrimitive, r[2], hZ, origin, rotation) -
        _get_N_prism(T,P,ClosedPrimitive, r[1], T(1.1) * hZ, origin, rotation) # increase volume to subtract
    else
        _get_N_prism(T,P,ClosedPrimitive, r, hZ, origin, rotation)
    end
    transform(g, transformations)
end

const PrismAliases = Dict{Int, String}(
    3 => "TriangularPrism", 
    4 => "QuadranglePrism", 
    5 => "PentagonalPrism", 
    6 => "HexagonalPrism"
)

function Dictionary(rp::RegularPrism{N,T, <:Any})::OrderedDict{String, Any} where {T, N}
    dict = OrderedDict{String, Any}()
    dict["r"] = rp.r # always a Real 
    dict["h"] = rp.hZ*2
    if !iszero(rp.origin) dict["origin"] = Dictionary(rp.origin) end
    if !isone(rp.rotation) dict["rotation"] = Dictionary(rp.rotation) end
    OrderedDict{String, Any}(PrismAliases[N] => dict)
end

function vertices(rp::RegularPrism{N,T,ClosedPrimitive,T}) where {T,N}
    xys = [rp.r .* sincos(T(2π)*(n-1)/N) for n in 1:N]
    pts = [CartesianPoint{T}(xy[2], xy[1], z) for z in (-rp.hZ, rp.hZ) for xy in xys]
    _transform_into_global_coordinate_system(pts, rp)
end
function vertices(rp::RegularPrism{N,T,OpenPrimitive,T}) where {T,N}
    xys = [rp.r .* sincos(T(2π)*(n-1)/N) for n in N:-1:1]
    pts = [CartesianPoint{T}(xy[2], xy[1], z) for z in (-rp.hZ, rp.hZ) for xy in xys]
    _transform_into_global_coordinate_system(pts, rp)
end

function surfaces(rp::RegularPrism{N,T,ClosedPrimitive,T}) where {T,N}
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


function sample(rp::RegularPrism{N,T,<:Any,T})::Vector{CartesianPoint{T}} where {T,N}
    [vertices(rp)...]
end

function _in(pt::CartesianPoint{T}, rp::RegularPrism{N,T,ClosedPrimitive,T}; csgtol::T = csg_default_tol(T)) where {T,N} 
    abs(pt.z) <= rp.hZ + csgtol && begin
        r, φ = hypot(pt.x, pt.y), atan(pt.y, pt.x)
        α = T(π/N)
        _r = r * cos(α - mod(φ, 2α)) / cos(α) 
        _r <= rp.r + csgtol
    end
end
function _in(pt::CartesianPoint{T}, rp::RegularPrism{N,T,OpenPrimitive,T}; csgtol::T = csg_default_tol(T)) where {T,N} 
    abs(pt.z) < rp.hZ - csgtol && begin
        r, φ = hypot(pt.x, pt.y), atan(pt.y, pt.x)
        α = T(π/N)
        _r = r * cos(α - mod(φ, 2α)) / cos(α) 
        _r < rp.r - csgtol
    end
end

extremum(rp::RegularPrism{N,T}) where {N,T} = hypot(rp.hZ, max((rp.r...)...))

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
