"""
    struct Torus{T,CO,TR,TP,TT,TT1,TT2} <: AbstractVolumePrimitive{T,CO}

Volume primitive describing a [Torus](@ref). It is defined as all points that are
within a given radius to a circle, parallel to the `xy` plane, with constant
radius around a given origin.

## Parametric types
* `T`: Precision type.
* `CO`: Describes whether the surface belongs to the primitive. 
    It can be `ClosedPrimitive`, i.e. the surface points belong to the primitive,
    or `OpenPrimitive`, i.e. the surface points do not belong to the primitive.
* `TR`: Type of `r_tube`.
    * `TR == T`: Full tube without cutout (constant radius `r_tube`).
* `TP`: Type of the azimuthial angle `φ`.
    * `TP == Nothing`: Full 2π in `φ`.
    * `TP == T`: Partial Torus ranging from `0` to `φ`.
* `TT`: Type of the polar angle `θ`.
    * `TT == Nothing`: Full 2π in `θ`.
    * `TP == Tuple{T,T}`: Partial Torus ranging from `θ[1]` to `θ[2]`.
* `TT1`: Type of the surface at `θ[1]` (`:inwards`, `:outwards` or `:flat`).
* `TT2`: Type of the surface at `θ[2]` (`:inwards`, `:outwards` or `:flat`).
    
## Fields
* `r_torus::T`: Distance of the center of the `Torus` to the center of the tube (in m).
* `r_tube::TR`: Radius of the tube of the `Torus` (in m).
* `φ::TP`: Range in azimuthial angle `φ` of the `Torus`.
* `θ::TT`: Range in polar angle `θ` of the `Torus`.
* `origin::CartesianPoint{T}`: The position of the center of the `Torus`.
* `rotation::SMatrix{3,3,T,9}`: Matrix that describes a rotation of the `Torus` around its `origin`.

## Definition in Configuration File

A `Torus` is defined in the configuration file as part of the `geometry` field 
of an object through the field `torus`.

Example definitions of a `Torus` looks like this:
```yaml
torus:
  r_torus: 10.0   # => r_torus = 10.0
  r_tube: 2       # => r_tube = 2.0
  phi: 
      from: 0.0°
      to: 360.0°  # => φ = nothing
  theta: 
      from: 0.0°
      to: 360.0°  # => θ = nothing
```
The fields `phi` and `theta` do not need to defined if they are full 2π.

Example definition of a `Torus` with an inner cut-out:
```yaml
torus:
  r_torus: 10.0   # => r_torus = 10.0
  r_tube: 
      from: 1.0
      to: 2.0     # => r_tube = (1.0, 2.0)
  phi: 
      from: 0.0°
      to: 360.0°  # => φ = nothing
  theta: 
      from: 0.0°
      to: 360.0°  # => θ = nothing
```
This is a `Torus` with `r_tube` having an inner radius of 1 and an outer radius of 2.

See also [Constructive Solid Geometry (CSG)](@ref).
"""
struct Torus{T,CO,TR,TP<:Union{Nothing,T},TT,TT1,TT2} <: AbstractVolumePrimitive{T,CO}
    r_torus::T
    r_tube::TR  # (r_tube_in, r_tube_out)
    φ::TP 
    θ::TT

    origin::CartesianPoint{T}
    rotation::SMatrix{3,3,T,9}
end

#Type conversion happens here
function Torus{T,CO}(r_torus, r_tube, φ::Tuple{T,T}, θ, origin, rotation) where {T,CO}
    r_torus = _csg_convert_args(T, r_torus)
    r_tube = _csg_convert_args(T, r_tube)
    φ = _csg_convert_args(T, φ)
    _φ = abs(φ[2]-φ[1])
    rotation=rotation*SMatrix{3}(cos(φ[1]),sin(φ[1]),0,-sin(φ[1]),cos(φ[1]),0,0,0,1)
    θ = _csg_convert_args(T, θ)
    TT1, TT2 = _get_conemantle_type(θ)
    Torus{T,CO,typeof(r_tube),typeof(_φ),typeof(θ),TT1,TT2}(r_torus, r_tube, _φ, θ, origin, rotation)
end

function Torus{T,CO}(r_torus, r_tube, φ, θ, origin, rotation) where {T,CO}
    r_torus = _csg_convert_args(T, r_torus)
    r_tube = _csg_convert_args(T, r_tube)
    φ = _csg_convert_args(T, φ)
    θ = _csg_convert_args(T, θ)
    TT1, TT2 = _get_conemantle_type(θ)
    Torus{T,CO,typeof(r_tube),typeof(φ),typeof(θ),TT1,TT2}(r_torus, r_tube, φ, θ, origin, rotation)
end

#Type promotion happens here
function Torus(CO, r_torus::TRTo, r_tube::TRTu, φ::TP, θ::TT, 
    origin::PT, rotation::ROT) where {TRTo, TRTu, TP, TT, PT, ROT}
        eltypes = _csg_get_promoted_eltype.((TRTo, TRTu, TP, TT, PT, ROT))
        T = float(promote_type(eltypes...))
        Torus{T,CO}(r_torus, r_tube, φ, θ, origin, rotation)
end

function Torus(::Type{CO} = ClosedPrimitive;
    r_torus = 1, 
    r_tube = 1,
    φ = nothing,
    θ = nothing,
    origin = zero(CartesianPoint{Int64}), 
    rotation = one(SMatrix{3, 3, Int64, 9})
) where {CO}
    Torus(CO, r_torus, r_tube, φ, θ, origin, rotation)
end

function Torus{T}(::Type{CO}=ClosedPrimitive;
    r_torus = 1., 
    r_tube = 1.,
    φ = nothing,
    θ = nothing,
    origin = zero(CartesianPoint{Float64}), 
    rotation = one(SMatrix{3, 3, Float64, 9})
) where {T, CO}
    Torus{T,CO}(r_torus, r_tube, φ, θ, origin, rotation)
end

Torus{T,CO,TR,TP,TT,TT1,TT2}( t::Torus{T,CO,TR,TP,TT}; COT = CO,
            origin::CartesianPoint{T} = t.origin,
            rotation::SMatrix{3,3,T,9} = t.rotation) where {T,CO<:Union{ClosedPrimitive, OpenPrimitive},TR,TP<:Union{Nothing,T},TT,TT1,TT2} =
    Torus{T,COT,TR,TP,TT,TT1,TT2}(t.r_torus, t.r_tube, t.φ, t.θ, origin, rotation)

const FullTorus{T,CO} = Torus{T,CO,T,Nothing,Nothing,Nothing,Nothing}
const FullPhiTorus{T,CO,TT1,TT2} = Torus{T,CO,T,Nothing,Tuple{T,T},TT1,TT2}
const FullThetaTorus{T,CO} = Torus{T,CO,T,T,Nothing,Nothing,Nothing}
const HollowTorus{T,CO} = Torus{T,CO,Tuple{T,T},Nothing,Nothing,Nothing,Nothing}
const HollowPhiTorus{T,CO,TT1,TT2} = Torus{T,CO,Tuple{T,T},Nothing,Tuple{T,T},TT1,TT2}
const HollowThetaTorus{T,CO} = Torus{T,CO,Tuple{T,T},T,Nothing,Nothing,Nothing}


function _get_conemantle_type(θ::Tuple{T,T})::Tuple{Symbol, Symbol} where {T}
    θ1::T, θ2::T = θ
    return ( mod(θ1, π) == 0 ? :flat : (mod(θ1, 2π) in 0..π ? :inwards : :outwards), 
             mod(θ2, π) == 0 ? :flat : (mod(θ2, 2π) in 0..π ? :outwards : :inwards))
end

_get_conemantle_type(::Nothing) = (Nothing, Nothing)


function Geometry(::Type{T}, ::Type{Torus}, dict::AbstractDict, input_units::NamedTuple, transformations::Transformations{T}) where {T}
    length_unit = input_units.length
    angle_unit = input_units.angle
    origin = get_origin(T, dict, length_unit)
    rotation = get_rotation(T, dict, angle_unit)

    r_torus = _parse_value(T, dict["r_torus"], length_unit)
    r_tube = _parse_radial_interval(T, dict["r_tube"], length_unit)
    r_tube isa Tuple{T,T} && iszero(r_tube[1]) ? r_tube = r_tube[2] : nothing
    φ_interval = parse_φ_of_primitive(T, dict, angle_unit)
    φ = if φ_interval isa Tuple{T,T}
        rotation = rotation * RotZ{T}(φ_interval[1])
        φ_interval[2] - φ_interval[1]
    elseif isnothing(φ_interval)
        nothing
    else
        throw(ConfigFileError("Error when trying to parse φ from configuration file."))
    end
    θ = parse_θ_of_primitive(T, dict, angle_unit)
    if haskey(dict, "z")
        @warn "Deprecation warning: Field `z` for `Torus` is deprecated. 
                Use instead (only) `origin` to specify the origin of the primitive.
                There might be a conflict with the possible field `origin`:
                The `z` component of the origin of the primitive is overwritten by the `z`."
        z = _parse_value(T, dict["z"], length_unit)
        origin = CartesianPoint{T}(origin[1], origin[2], z)
    end

    t = Torus{T}(ClosedPrimitive;
            r_torus = r_torus, r_tube = r_tube, φ = φ, θ = θ, 
            origin = origin, rotation = rotation)   
    transform(t, transformations)
end

function Dictionary(t::Torus{T,<:Any,TR})::OrderedDict{String, Any} where {T,TR}
    dict = OrderedDict{String, Any}()
    dict["r_torus"] = t.r_torus
    dict["r_tube"] = TR <: Real ? t.r_tube : OrderedDict{String, Any}("from" => t.r_tube[1], "to" => t.r_tube[2]) 
    if !isnothing(t.φ) dict["phi"]   = OrderedDict("from" => "0°", "to" => string(rad2deg(t.φ))*"°") end
    if !isnothing(t.θ) dict["theta"] = OrderedDict("from" => string(rad2deg(t.θ[1]))*"°", "to" => string(rad2deg(t.θ[2]))*"°") end
    if t.origin != zero(CartesianVector{T}) dict["origin"] = t.origin end
    if t.rotation != one(SMatrix{3,3,T,9}) 
        d = Dictionary(t.rotation)
        if unique(keys(d)) == ["Z"]
            φ0 = mod2pi(_parse_value(T, d["Z"], internal_angle_unit))
            dict["phi"] = OrderedDict("from" => string(rad2deg(φ0))*"°", "to" => string(rad2deg(φ0 + t.φ))*"°")
        else
            dict["rotation"] = d
        end
    end
    OrderedDict{String, Any}("torus" => dict)
end

function surfaces(t::FullTorus{T,ClosedPrimitive}) where {T}
    tm = FullTorusMantle{T,:inwards}(t.r_torus, t.r_tube, t.φ, t.θ, t.origin, t.rotation)
    (tm, )
end
function surfaces(t::FullTorus{T,OpenPrimitive}) where {T}
    tm = FullTorusMantle{T,:outwards}(t.r_torus, t.r_tube, t.φ, t.θ, t.origin, t.rotation)
    (tm, )
end
function surfaces(t::FullPhiTorus{T,ClosedPrimitive,TT1,TT2}) where {T,TT1,TT2}
    tm = FullPhiTorusMantle{T,:inwards}(t.r_torus, t.r_tube, t.φ, t.θ, t.origin, t.rotation)
    θ1, θ2 = t.θ
    cm1 = TorusThetaSurface((t.r_torus, t.r_torus + t.r_tube * cos(θ1)), t.φ, t.r_tube * sin(θ1)/2, t.origin + t.rotation * CartesianVector{T}(0,0,t.r_tube * sin(θ1)/2), t.rotation, Val{TT1}())
    cm2 = TorusThetaSurface((t.r_torus, t.r_torus + t.r_tube * cos(θ2)), t.φ, t.r_tube * sin(θ2)/2, t.origin + t.rotation * CartesianVector{T}(0,0,t.r_tube * sin(θ2)/2), t.rotation, Val{TT2}())
    (tm, cm1, cm2)
end
function surfaces(t::FullPhiTorus{T,OpenPrimitive,TT1,TT2}) where {T,TT1,TT2}
    tm = FullPhiTorusMantle{T,:outwards}(t.r_torus, t.r_tube, t.φ, t.θ, t.origin, t.rotation)
    θ1, θ2 = t.θ
    cm1 = flip(TorusThetaSurface((t.r_torus, t.r_torus + t.r_tube * cos(θ1)), t.φ, t.r_tube * sin(θ1)/2, t.origin + t.rotation * CartesianVector{T}(0,0,t.r_tube * sin(θ1)/2), t.rotation, Val{TT1}()))
    cm2 = flip(TorusThetaSurface((t.r_torus, t.r_torus + t.r_tube * cos(θ2)), t.φ, t.r_tube * sin(θ2)/2, t.origin + t.rotation * CartesianVector{T}(0,0,t.r_tube * sin(θ2)/2), t.rotation, Val{TT2}()))
    (tm, cm1, cm2)
end
function surfaces(t::FullThetaTorus{T,ClosedPrimitive}) where {T}
    tm = FullThetaTorusMantle{T,:inwards}(t.r_torus, t.r_tube, t.φ, t.θ, t.origin, t.rotation)
    es1 = CircularArea{T}(t.r_tube, t.θ, t.origin + t.rotation * CartesianVector{T}(t.r_torus, 0, 0), t.rotation * RotX(-π/2) )
    es2 = CircularArea{T}(t.r_tube, t.θ, t.origin + t.rotation * CartesianVector{T}(t.r_torus * cos(t.φ), t.r_torus * sin(t.φ), 0), t.rotation * RotZ(t.φ) * RotX(π/2) )
    (tm, es1, es2)
end
function surfaces(t::FullThetaTorus{T,OpenPrimitive}) where {T}
    tm = FullThetaTorusMantle{T,:outwards}(t.r_torus, t.r_tube, t.φ, t.θ, t.origin, t.rotation)
    φ1, φ2 = T(0), t.φ
    es1 = CircularArea{T}(t.r_tube, t.θ, t.origin + t.rotation * CartesianVector{T}(t.r_torus, 0, 0), t.rotation * RotX(π/2) )
    es2 = CircularArea{T}(t.r_tube, t.θ, t.origin + t.rotation * CartesianVector{T}(t.r_torus * cos(t.φ), t.r_torus * sin(t.φ), 0), t.rotation * RotZ(t.φ) * RotX(-π/2) )
    (tm, es1, es2)
end
function surfaces(t::Torus{T,ClosedPrimitive,T,T,Tuple{T,T},TT1,TT2}) where {T,TT1,TT2}
    tm = TorusMantle{T,T,Tuple{T,T},:inwards}(t.r_torus, t.r_tube, t.φ, t.θ, t.origin, t.rotation)
    θ1, θ2 = t.θ
    cm1 = TorusThetaSurface((t.r_torus, t.r_torus + t.r_tube * cos(θ1)), t.φ, t.r_tube * sin(θ1)/2, t.origin + t.rotation * CartesianVector{T}(0,0,t.r_tube * sin(θ1)/2), t.rotation, Val{TT1}())
    cm2 = TorusThetaSurface((t.r_torus, t.r_torus + t.r_tube * cos(θ2)), t.φ, t.r_tube * sin(θ2)/2, t.origin + t.rotation * CartesianVector{T}(0,0,t.r_tube * sin(θ2)/2), t.rotation, Val{TT2}())
    es1 = flip(PartialCircularArea{T}(t.r_tube, abs(-(t.θ...)), t.origin + t.rotation * CartesianVector{T}(t.r_torus, 0, 0), t.rotation * RotX(π/2) * RotZ{T}(t.θ[1])))
    es2 = PartialCircularArea{T}(t.r_tube, abs(-(t.θ...)), t.origin + t.rotation * CartesianVector{T}(t.r_torus * cos(t.φ), t.r_torus * sin(t.φ), 0), t.rotation * RotZ(t.φ) * RotX(π/2) * RotZ{T}(t.θ[1]))
    (tm, cm1, cm2, es1, es2)
end
function surfaces(t::Torus{T,OpenPrimitive,T,T,Tuple{T,T},TT1,TT2}) where {T,TT1,TT2}
    tm = TorusMantle{T,T,Tuple{T,T},:outwards}(t.r_torus, t.r_tube, t.φ, t.θ, t.origin, t.rotation)
    θ1, θ2 = t.θ
    cm1 = flip(TorusThetaSurface((t.r_torus, t.r_torus + t.r_tube * cos(θ1)), t.φ, t.r_tube * sin(θ1)/2, t.origin + t.rotation * CartesianVector{T}(0,0,t.r_tube * sin(θ1)/2), t.rotation, Val{TT1}()))
    cm2 = flip(TorusThetaSurface((t.r_torus, t.r_torus + t.r_tube * cos(θ2)), t.φ, t.r_tube * sin(θ2)/2, t.origin + t.rotation * CartesianVector{T}(0,0,t.r_tube * sin(θ2)/2), t.rotation, Val{TT2}()))
    es1 = PartialCircularArea{T}(t.r_tube, abs(-(t.θ...)), t.origin + t.rotation * CartesianVector{T}(t.r_torus, 0, 0), t.rotation * RotX(π/2) * RotZ{T}(t.θ[1]))
    es2 = flip(PartialCircularArea{T}(t.r_tube, abs(-(t.θ...)), t.origin + t.rotation * CartesianVector{T}(t.r_torus * cos(t.φ), t.r_torus * sin(t.φ), 0), t.rotation * RotZ(t.φ) * RotX(π/2) * RotZ{T}(t.θ[1])))
    (tm, cm1, cm2, es1, es2)
end

function surfaces(t::HollowTorus{T,ClosedPrimitive}) where {T}
    tm_in = FullTorusMantle{T,:outwards}(t.r_torus, t.r_tube[1], t.φ, t.θ, t.origin, t.rotation)
    tm_out = FullTorusMantle{T,:inwards}(t.r_torus, t.r_tube[2], t.φ, t.θ, t.origin, t.rotation)
    (tm_in, tm_out)
end
function surfaces(t::HollowTorus{T,OpenPrimitive}) where {T}
    tm_in = FullTorusMantle{T,:inwards}(t.r_torus, t.r_tube[1], t.φ, t.θ, t.origin, t.rotation)
    tm_out = FullTorusMantle{T,:outwards}(t.r_torus, t.r_tube[2], t.φ, t.θ, t.origin, t.rotation)
    (tm_in, tm_out)
end
function surfaces(t::HollowPhiTorus{T,ClosedPrimitive,TT1,TT2}) where {T,TT1,TT2}
    tm_in = FullPhiTorusMantle{T,:outwards}(t.r_torus, t.r_tube[1], t.φ, t.θ, t.origin, t.rotation)
    tm_out = FullPhiTorusMantle{T,:inwards}(t.r_torus, t.r_tube[2], t.φ, t.θ, t.origin, t.rotation)
    θ1, θ2 = t.θ
    cm1 = TorusThetaSurface(t.r_torus .+ t.r_tube .* cos(θ1), t.φ, abs(-(t.r_tube...)) * sin(θ1)/2, t.origin + t.rotation * CartesianVector{T}(0,0,mean(t.r_tube) * sin(θ1)), t.rotation, Val{TT1}())
    cm2 = TorusThetaSurface(t.r_torus .+ t.r_tube .* cos(θ2), t.φ, abs(-(t.r_tube...)) * sin(θ2)/2, t.origin + t.rotation * CartesianVector{T}(0,0,mean(t.r_tube) * sin(θ2)), t.rotation, Val{TT2}())
    (tm_in, tm_out, cm1, cm2)
end
function surfaces(t::HollowPhiTorus{T,OpenPrimitive,TT1,TT2}) where {T,TT1,TT2}
    tm_in = FullPhiTorusMantle{T,:inwards}(t.r_torus, t.r_tube[1], t.φ, t.θ, t.origin, t.rotation)
    tm_out = FullPhiTorusMantle{T,:outwards}(t.r_torus, t.r_tube[2], t.φ, t.θ, t.origin, t.rotation)
    θ1, θ2 = t.θ
    cm1 = flip(TorusThetaSurface(t.r_torus .+ t.r_tube .* cos(θ1), t.φ, abs(-(t.r_tube...)) * sin(θ1)/2, t.origin + t.rotation * CartesianVector{T}(0,0,mean(t.r_tube) * sin(θ1)), t.rotation, Val{TT1}()))
    cm2 = flip(TorusThetaSurface(t.r_torus .+ t.r_tube .* cos(θ2), t.φ, abs(-(t.r_tube...)) * sin(θ2)/2, t.origin + t.rotation * CartesianVector{T}(0,0,mean(t.r_tube) * sin(θ2)), t.rotation, Val{TT2}()))
    (tm_in, tm_out, cm1, cm2)
end
function surfaces(t::HollowThetaTorus{T,ClosedPrimitive}) where {T}
    tm_in = FullThetaTorusMantle{T,:outwards}(t.r_torus, t.r_tube[1], t.φ, t.θ, t.origin, t.rotation)
    tm_out = FullThetaTorusMantle{T,:inwards}(t.r_torus, t.r_tube[2], t.φ, t.θ, t.origin, t.rotation)
    es1 = Annulus{T}(t.r_tube, abs(-(t.θ...)), t.origin + t.rotation * CartesianVector{T}(t.r_torus, 0, 0), t.rotation * RotX(-π/2)  * RotZ{T}(t.θ[1]))
    es2 = Annulus{T}(t.r_tube, abs(-(t.θ...)), t.origin + t.rotation * CartesianVector{T}(t.r_torus * cos(t.φ), t.r_torus * sin(t.φ), 0), t.rotation * RotZ(t.φ) * RotX(π/2)  * RotZ{T}(t.θ[1]))
    (tm_in, tm_out, es1, es2)
end
function surfaces(t::HollowThetaTorus{T,OpenPrimitive}) where {T}
    tm_in = FullThetaTorusMantle{T,:inwards}(t.r_torus, t.r_tube[1], t.φ, t.θ, t.origin, t.rotation)
    tm_out = FullThetaTorusMantle{T,:outwards}(t.r_torus, t.r_tube[2], t.φ, t.θ, t.origin, t.rotation)
    es1 = Annulus{T}(t.r_tube, abs(-(t.θ...)), t.origin + t.rotation * CartesianVector{T}(t.r_torus, 0, 0), t.rotation * RotX(π/2)  * RotZ{T}(t.θ[1]))
    es2 = Annulus{T}(t.r_tube, abs(-(t.θ...)), t.origin + t.rotation * CartesianVector{T}(t.r_torus * cos(t.φ), t.r_torus * sin(t.φ), 0), t.rotation * RotZ(t.φ) * RotX(-π/2)  * RotZ{T}(t.θ[1]))
    (tm_in, tm_out, es1, es2)
end
function surfaces(t::Torus{T,ClosedPrimitive,Tuple{T,T},T,Tuple{T,T},TT1,TT2}) where {T,TT1,TT2}
    tm_in = TorusMantle{T,T,Tuple{T,T},:outwards}(t.r_torus, t.r_tube[1], t.φ, t.θ, t.origin, t.rotation)
    tm_out = TorusMantle{T,T,Tuple{T,T},:inwards}(t.r_torus, t.r_tube[2], t.φ, t.θ, t.origin, t.rotation)
    θ1, θ2 = t.θ
    cm1 = TorusThetaSurface(t.r_torus .+ t.r_tube .* cos(θ1), t.φ, abs(-(t.r_tube...)) * sin(θ1)/2, t.origin + t.rotation * CartesianVector{T}(0,0,mean(t.r_tube) * sin(θ1)), t.rotation, Val{TT1}())
    cm2 = TorusThetaSurface(t.r_torus .+ t.r_tube .* cos(θ2), t.φ, abs(-(t.r_tube...)) * sin(θ2)/2, t.origin + t.rotation * CartesianVector{T}(0,0,mean(t.r_tube) * sin(θ2)), t.rotation, Val{TT2}())
    es1 = flip(PartialAnnulus{T}(t.r_tube, abs(-(t.θ...)), t.origin + t.rotation * CartesianVector{T}(t.r_torus, 0, 0), t.rotation * RotX(π/2) * RotZ{T}(t.θ[1])))
    es2 = PartialAnnulus{T}(t.r_tube, abs(-(t.θ...)), t.origin + t.rotation * CartesianVector{T}(t.r_torus * cos(t.φ), t.r_torus * sin(t.φ), 0), t.rotation * RotZ(t.φ) * RotX(π/2) * RotZ{T}(t.θ[1]))
    (tm_in, tm_out, cm1, cm2, es1, es2)
end
function surfaces(t::Torus{T,OpenPrimitive,Tuple{T,T},T,Tuple{T,T},TT1,TT2}) where {T,TT1,TT2}
    tm_in = TorusMantle{T,T,Tuple{T,T},:inwards}(t.r_torus, t.r_tube[1], t.φ, t.θ, t.origin, t.rotation)
    tm_out = TorusMantle{T,T,Tuple{T,T},:outwards}(t.r_torus, t.r_tube[2], t.φ, t.θ, t.origin, t.rotation)
    θ1, θ2 = t.θ
    cm1 = flip(TorusThetaSurface(t.r_torus .+ t.r_tube .* cos(θ1), t.φ, abs(-(t.r_tube...)) * sin(θ1)/2, t.origin + t.rotation * CartesianVector{T}(0,0,mean(t.r_tube) * sin(θ1)), t.rotation, Val{TT1}()))
    cm2 = flip(TorusThetaSurface(t.r_torus .+ t.r_tube .* cos(θ2), t.φ, abs(-(t.r_tube...)) * sin(θ2)/2, t.origin + t.rotation * CartesianVector{T}(0,0,mean(t.r_tube) * sin(θ2)), t.rotation, Val{TT2}()))
    es1 = PartialAnnulus{T}(t.r_tube, abs(-(t.θ...)), t.origin + t.rotation * CartesianVector{T}(t.r_torus, 0, 0), t.rotation * RotX(π/2) * RotZ{T}(t.θ[1]))
    es2 = flip(PartialAnnulus{T}(t.r_tube, abs(-(t.θ...)), t.origin + t.rotation * CartesianVector{T}(t.r_torus * cos(t.φ), t.r_torus * sin(t.φ), 0), t.rotation * RotZ(t.φ) * RotX(π/2) * RotZ{T}(t.θ[1])))
    (tm_in, tm_out, cm1, cm2, es1, es2)
end


function _in(pt::CartesianPoint{T}, t::Torus{T,ClosedPrimitive}; csgtol::T = csg_default_tol(T)) where {T}
    _r = hypot(hypot(pt.x, pt.y) - t.r_torus, pt.z)
    rmin::T, rmax::T = _radial_endpoints(t.r_tube)
    return rmin - csgtol <= _r <= rmax + csgtol &&
        (isnothing(t.φ) || _in_angular_interval_closed(atan(pt.y, pt.x), t.φ, csgtol = csgtol)) &&
        (isnothing(t.θ) || _in_angular_interval_closed(atan(pt.z, hypot(pt.x, pt.y) - t.r_torus), t.θ, csgtol = csgtol))
end
function _in(pt::CartesianPoint{T}, t::Torus{T,OpenPrimitive}; csgtol::T = csg_default_tol(T)) where {T}
    _r = hypot(hypot(pt.x, pt.y) - t.r_torus, pt.z)
    rmin::T, rmax::T = _radial_endpoints(t.r_tube)
    return rmin + csgtol < _r < rmax - csgtol &&
        (isnothing(t.φ) || _in_angular_interval_open(atan(pt.y, pt.x), t.φ, csgtol = csgtol)) &&
        (isnothing(t.θ) || _in_angular_interval_open(atan(pt.z, hypot(pt.x, pt.y) - t.r_torus), t.θ, csgtol = csgtol))
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


extremum(t::Torus{T}) where {T} = t.r_torus + max(t.r_tube...)

# function Dictionary(t::Torus{T}) where {T}
#     dict = OrderedDict{String,Any}()
#     dict["r_torus"] = t.r_torus
#     dict["r_tube"] = typeof(t.r_tube) == T ? t.r_tube : OrderedDict{String,Any}("from" => t.r_tube.left, "to" => t.r_tube.right)
#     if !isnothing(t.φ) dict["phi"] = OrderedDict{String,Any}("from" => t.φ.left, "to" => t.φ.right) end
#     if !isnothing(t.θ) dict["theta"] = OrderedDict{String,Any}("from" => t.θ.left, "to" => t.θ.right) end
#     if t.z != 0 dict["z"] = t.z end
#     OrderedDict{String,Any}("torus" => dict)
# end


# get_r_tube_limits(t::Torus{T}) where {T} = _radial_endpoints(t.r_tube)
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
