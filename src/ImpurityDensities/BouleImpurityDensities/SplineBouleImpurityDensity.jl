"""
    struct LinBouleImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

a + b*z
Impurity density model which assumes a linear gradient in impurity density in z.
 
## Fields
* `a:T`: impurity density values at the boule origin.
* `b::T`: linear slope in `z` direction.
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
    step = z[2] - z[1]
    ρ = SVector{N, T}(to_internal_units.(ρ))
    spline = cubic_spline_interpolation(range(z[1], stop = z[end], length = N), ρ)
    SplineBouleImpurityDensity{N, T}(SVector{N, T}(z), ρ, spline, T(to_internal_units(det_z0)))
end

function get_impurity_density(idm::SplineBouleImpurityDensity, pt::AbstractCoordinatePoint{T})::T where {T}
    cpt = CartesianPoint(pt)
    z = cpt[3]
    idm.spline(idm.det_z0 - z)
end

(*)(scale::Real, idm::SplineBouleImpurityDensity{N, T}) where {N, T} = SplineBouleImpurityDensity{T}(Vector(idm.z), Vector(scale*idm.ρ), idm.det_z0)

(+)(offset::Union{<:Real, <:Quantity{<:Real, Unitful.𝐋^(-3)}}, idm::SplineBouleImpurityDensity{N, T}) where {N, T} = SplineBouleImpurityDensity{T}(Vector(idm.z), Vector(idm.ρ .+ to_internal_units(offset)), idm.det_z0)