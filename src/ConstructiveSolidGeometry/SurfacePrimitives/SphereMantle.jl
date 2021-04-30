struct SphereMantle{T} <: AbstractSurfacePrimitive{T}
    r::T
end

in(p::AbstractCoordinatePoint, s::SphereMantle) = _isapprox_sph_r(p, s.r)

#=
function sample(s::SphereMantle{T}, step::Real)::Vector{CylindricalPoint{T}} where {T}
    samples = [
        CylindricalPoint{T}(sqrt(s.r^2 - z^2),φ,z)
        for z in in -s.r:step:s.r
        for φ in 0:step/s.r:2π
    ]
end
=#

function sample(s::SphereMantle{T}, Nsamps::NTuple{3,Int})::Vector{CylindricalPoint{T}} where {T}
    samples = [
        CylindricalPoint{T}(sqrt(s.r^2 - z^2),φ,z)
        for z in (Nsamps[3] ≤ 1 ? s.r : range(-s.r, s.r, length = Nsamps[3]))
        for φ in (Nsamps[2] ≤ 1 ? 0 : range(0, 2π, length = Nsamps[2]))
    ]
end
