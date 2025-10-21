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

function SplineBouleImpurityDensity{T}(z::Vector{<:Number}, ρ::Vector{<:Number}, det_z0::Number) where {T}
    N = length(z)
    @assert N == length(ρ) > 1 "Vectors must be the same length, and at least 2 points are required for cubic spline interpolation."
    z = to_internal_units.(z)
    ρ = SVector{N, T}(to_internal_units.(ρ))
    spline = cubic_spline_interpolation(range(z[1], stop = z[end], length = N), ρ)
    SplineBouleImpurityDensity{N, T}(SVector{N, T}(z), ρ, spline, T(to_internal_units(det_z0)))
end

function get_impurity_density(idm::SplineBouleImpurityDensity, pt::AbstractCoordinatePoint{T})::T where {T}
    cpt = CartesianPoint(pt)
    z = cpt[3]
    idm.spline(idm.det_z0 - z)
end

function resample_with_min_step(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    @assert length(x) == length(y) "x and y must have the same length"
    @assert issorted(x) "x must be sorted in increasing order"

    # Compute the minimum step
    dxs = diff(x)
    min_step = minimum(dxs)

    # Construct StepRangeLen from min step
    x_min, x_max = first(x), last(x)
    n = Int(cld(x_max - x_min, min_step)) + 1  # cld = ceil division
    new_x = StepRangeLen(x_min, min_step, n)

    # Interpolate with flat extrapolation
    itp = LinearInterpolation(x, y, extrapolation_bc=Flat())
    new_y = SVector{n}(itp.(new_x))

    new_x, new_y
end

(*)(scale::Real, idm::SplineBouleImpurityDensity{N, T}) where {N, T} = SplineBouleImpurityDensity{T}(Vector(idm.z), Vector(scale*idm.ρ), idm.det_z0)

(+)(offset::Union{<:Real, <:Quantity{<:Real, Unitful.𝐋^(-3)}}, idm::SplineBouleImpurityDensity{N, T}) where {N, T} = SplineBouleImpurityDensity{T}(Vector(idm.z), Vector(idm.ρ .+ to_internal_units(offset)), idm.det_z0)