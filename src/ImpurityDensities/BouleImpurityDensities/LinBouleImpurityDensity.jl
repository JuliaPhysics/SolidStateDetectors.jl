"""
    struct LinBouleImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

a + b*z
Impurity density model which assumes a linear gradient in impurity density in z.
 
## Fields
* `a:T`: impurity density values at the boule origin.
* `b::T`: linear slope in `z` direction.
* `det_z0::T`: z coordinate of the detector origin in boule coordinates. The z-direction of the detector is opposite to the z-direction of the boule coordinates.
"""

struct LinBouleImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}
    a::T
    b::T 
    det_z0::T
end

function LinBouleImpurityDensity{T}(pars::Vector{<:Number}, det_z0::Number) where {T}
    LinBouleImpurityDensity{T}(
        T.(to_internal_units(pars))..., 
        T(to_internal_units(det_z0))
        )
end

function ImpurityDensity(T::DataType, t::Val{:linear_boule}, dict::AbstractDict, input_units::NamedTuple)
    a::T = haskey(dict, "a") ? _parse_value(T, dict["a"], input_units.length^(-3)) : T(0)
    b::T = haskey(dict, "b") ? _parse_value(T, dict["b"], input_units.length^(-4)) : T(0)
    det_z0::T = haskey(dict, "det_z0") ? _parse_value(T, dict["det_z0"], input_units.length) : T(0)
    LinBouleImpurityDensity{T}(a, b, det_z0)
end

function get_impurity_density(idm::LinBouleImpurityDensity, pt::AbstractCoordinatePoint{T})::T where {T}
    cpt = CartesianPoint(pt)
    z = cpt[3]
    idm.a + idm.b * (idm.det_z0 - z)
end

(*)(scale::Real, idm::LinBouleImpurityDensity{T}) where {T} = LinBouleImpurityDensity{T}(T(scale*idm.a), T(scale*idm.b), idm.det_z0)

(+)(offset::Union{<:Real, <:Quantity{<:Real, Unitful.𝐋^(-3)}}, idm::LinBouleImpurityDensity{T}) where {T} = LinBouleImpurityDensity{T}(T(to_internal_units(offset)+idm.a), idm.b, idm.det_z0)