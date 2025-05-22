"""
    struct ParExpBouleImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

a + b*z + c*z^2 + n*exp((z-l)/m)
Impurity density model which assumes a linear + exponential component in impurity density in z.

## Fields
* `a:T`: impurity density values at the boule origin.
* `b::T`: linear slope in `z` direction.
* `c::T`: quadratic coefficient.
* `n::T`: exponential coefficient.
* `l::T`: exponential offset.
* `m::T`: exponential scale factor.
* `det_z0::T`: z coordinate of the detector origin in boule coordinates. The z-direction of the detector is opposite to the z-direction of the boule coordinates.
"""

struct ParExpBouleImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}
    a::T
    b::T 
    c::T 
    n::T
    l::T
    m::T
    det_z0::T
end

function ParExpBouleImpurityDensity{T}(pars::Vector{Number}, det_z0::Number) where {T}
    ParExpBouleImpurityDensity{T}(
        T.(to_internal_units(pars))..., 
        T(to_internal_units(det_z0))
        )
end

function ImpurityDensity(T::DataType, t::Val{:parabolic_exponential_boule}, dict::AbstractDict, input_units::NamedTuple)
    a::T = haskey(dict, "a") ? _parse_value(T, dict["a"], input_units.length^(-3)) : T(0)
    b::T = haskey(dict, "b") ? _parse_value(T, dict["b"], input_units.length^(-4)) : T(0)
    c::T = haskey(dict, "c") ? _parse_value(T, dict["c"], input_units.length^(-5)) : T(0)
    n::T = haskey(dict, "n") ? _parse_value(T, dict["n"], input_units.length^(-3)) : T(0)
    l::T = haskey(dict, "l") ? _parse_value(T, dict["l"], input_units.length) : T(0)
    m::T = haskey(dict, "m") ? _parse_value(T, dict["m"], input_units.length) : T(1)
    det_z0::T = haskey(dict, "det_z0") ? _parse_value(T, dict["det_z0"], input_units.length) : T(0)
    ParExpBouleImpurityDensity{T}(a, b, c, n, l, m, det_z0)
end

function get_impurity_density(idm::ParExpBouleImpurityDensity, pt::AbstractCoordinatePoint{T})::T where {T}
    cpt = CartesianPoint(pt)
    z = cpt[3]
    idm.a + idm.b * (idm.det_z0 - z) + idm.c * (idm.det_z0 - z)^2 + idm.n * exp((idm.det_z0 - z - idm.l)/idm.m)
end

(*)(scale::Real, idm::ParExpBouleImpurityDensity{T}) where {T} = ParExpBouleImpurityDensity{T}(T(scale*idm.a), T(scale*idm.b), T(scale*idm.c), T(scale*idm.n), idm.l, idm.m, idm.det_z0)

(+)(offset::Union{<:Real, <:Quantity{<:Real, Unitful.𝐋^(-3)}}, idm::ParExpBouleImpurityDensity{T}) where {T} = ParExpBouleImpurityDensity{T}(T(to_internal_units(offset)+idm.a), idm.b, idm.c, idm.n, idm.l, idm.m, idm.det_z0)