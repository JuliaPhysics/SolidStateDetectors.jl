"""
    struct SplineBouleImpurityDensity{N, T <: SSDFloat} <: AbstractImpurityDensity{T}

Impurity density model which assumes an arbitrary impurity density in z (boule coordinates) defined by cubic spline interpolation between given points.

## Fields
* `z::SVector{N, T}`: z coordinates (in boule coordinates) where the impurity density is defined.
* `ρ::SVector{N, T}`: impurity density values at the given z coordinates.
* `spline::Interpolations.Extrapolation{T}`: cubic spline interpolation of the impurity density values.
* `det_z0::T`: z coordinate of the detector origin in boule coordinates. The z-direction of the detector is opposite to the z-direction of the boule coordinates.
"""
struct SplineBouleImpurityDensity{N, T <: SSDFloat} <: AbstractImpurityDensity{T}
    z::SVector{N, T}
    ρ::SVector{N, T}
    spline::Interpolations.Extrapolation{T}
    det_z0::T
end

function SplineBouleImpurityDensity{T}(z::AbstractVector{<:RealQuantity}, ρ::AbstractVector{<:RealQuantity}, det_z0::RealQuantity) where {T}
    N = length(z)
    @assert N == length(ρ) > 1 "Vectors must be the same length, and at least 2 points are required for cubic spline interpolation."
    
    idx = sortperm(z)
    z = SVector{N, T}(to_internal_units.(z[idx]))
    ρ = SVector{N, T}(to_internal_units.(ρ[idx]))

    dxs = diff(z)
    z_range, resampled_ρ = if all(x -> isapprox(x, dxs[1]), dxs)
        range(first(z), stop = last(z), length = length(z)), ρ
    else
        @info "Resampling impurity density to uniform steps in z"
        min_step = minimum(dxs)
        z_min, z_max = first(z), last(z)
        steps = (z_max - z_min)/min_step
        n = isapprox(ceil(steps) - steps, 1) ? Int(ceil(steps)) : Int(ceil(steps)) + 1
        z_r = range(z_min, stop = z_max, length = n)
        itp = Interpolations.LinearInterpolation(z, ρ, extrapolation_bc=Flat())
        @info steps, n, z_r
        z_r, itp.(z_r)
    end

    spline = Interpolations.cubic_spline_interpolation(z_range, resampled_ρ)
    SplineBouleImpurityDensity{N, T}(z, ρ, spline, T(to_internal_units(det_z0)))
end

function get_impurity_density(idm::SplineBouleImpurityDensity, pt::AbstractCoordinatePoint{T})::T where {T}
    idm.spline(idm.det_z0 - pt.z)
end

(*)(scale::Real, idm::SplineBouleImpurityDensity{N, T}) where {N, T} = SplineBouleImpurityDensity{T}(idm.z, scale*idm.ρ, idm.det_z0)

(+)(offset::Union{<:Real, <:Quantity{<:Real, Unitful.𝐋^(-3)}}, idm::SplineBouleImpurityDensity{N, T}) where {N, T} = SplineBouleImpurityDensity{T}(idm.z, idm.ρ .+ to_internal_units(offset), idm.det_z0)