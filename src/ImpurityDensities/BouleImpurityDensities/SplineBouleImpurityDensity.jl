"""
    struct SplineBouleImpurityDensity{N, T <: SSDFloat} <: AbstractImpurityDensity{T}

Impurity density model which assumes an arbitrary impurity density in z (boule coordinates) defined by FritschCarlsonMonotonicInterpolation between given points.

## Fields
* `z::SVector{N, T}`: z coordinates (in boule coordinates) where the impurity density is defined.
* `ρ::SVector{N, T}`: impurity density values at the given z coordinates.
* `spline::Interpolations.MonotonicInterpolation{T}`: FritschCarlsonMonotonicInterpolation of the impurity density values.
* `det_z0::T`: z coordinate of the detector origin in boule coordinates. The z-direction of the detector is opposite to the z-direction of the boule coordinates.
"""
struct SplineBouleImpurityDensity{N, T <: SSDFloat} <: AbstractImpurityDensity{T}
    z::SVector{N, T}
    ρ::SVector{N, T}
    spline::Interpolations.MonotonicInterpolation{T}
    det_z0::T
end

function SplineBouleImpurityDensity{T}(z::AbstractVector{<:RealQuantity}, ρ::AbstractVector{<:RealQuantity}, det_z0::RealQuantity) where {T}
    N = length(z)
    @assert N == length(ρ) > 1 "Vectors must be the same length, and at least 2 points are required for cubic spline interpolation."
    
    idx = sortperm(z)
    z = SVector{N, T}(to_internal_units.(z[idx]))
    ρ = SVector{N, T}(to_internal_units.(ρ[idx]))

    spline = interpolate(z, ρ, FritschCarlsonMonotonicInterpolation())
    SplineBouleImpurityDensity{N, T}(z, ρ, spline, T(to_internal_units(det_z0)))
end

function get_impurity_density(idm::SplineBouleImpurityDensity, pt::AbstractCoordinatePoint{T})::T where {T}
    idm.spline(idm.det_z0 - pt.z)
end

(*)(scale::Real, idm::SplineBouleImpurityDensity{N, T}) where {N, T} = SplineBouleImpurityDensity{T}(idm.z, scale*idm.ρ, idm.det_z0)

(+)(offset::Union{<:Real, <:Quantity{<:Real, Unitful.𝐋^(-3)}}, idm::SplineBouleImpurityDensity{N, T}) where {N, T} = SplineBouleImpurityDensity{T}(idm.z, idm.ρ .+ to_internal_units(offset), idm.det_z0)