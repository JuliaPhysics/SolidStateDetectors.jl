struct Torus{T,TR,TB,TP,TT} <: AbstractVolumePrimitive{T}
    r_torus::TR
    r_tube::TB
    φ::TP
    θ::TT
    function Torus( ::Type{T},
                   r_torus::T,
                   r_tube::Union{T, <:AbstractInterval{T}},
                   φ::Union{Nothing, <:AbstractInterval{T}},
                   θ::Union{Nothing, <:AbstractInterval{T}}) where {T}
        new{T,typeof(r_torus),typeof(r_tube),typeof(φ),typeof(θ)}(r_torus, r_tube, φ, θ)
    end
end

#Constructors
function Torus(;r_torus = 1, r_tubeMin = 0, r_tubeMax = 1, φMin = 0, φMax = 2π, θMin = 0, θMax = 2π)
    T = float(promote_type(typeof.((r_torus, r_tubeMin, r_tubeMax, φMin, φMax, θMin, θMax))...))
    r_torus = T(r_torus)
    r_tube = r_tubeMin == 0 ? T(r_tubeMax) : T(r_tubeMin)..T(r_tubeMax)
    φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
    θ = mod(T(θMax) - T(θMin), T(2π)) == 0 ? nothing : T(θMin)..T(θMax)
    Torus( T, r_torus, r_tube, φ, θ)
end
Torus(r_torus, r_tubeMin, r_tubeMax, φMin, φMax, θMin, θMax) = Torus(;r_torus, r_tubeMin, r_tubeMax, φMin, φMax, θMin, θMax)


_in_torr_r_tube(p::CartesianPoint, r_torus::Real, r_tube::Real) = hypot(hypot(p.x, p.y) - r_torus, p.z) <= r_tube
_in_torr_r_tube(p::CartesianPoint, r_torus::Real, r_tube::AbstractInterval) = hypot(hypot(p.x, p.y) - r_torus, p.z) in r_tube

_in_torr_θ(p::CartesianPoint{T}, r_torus::Real, θ::AbstractInterval) where {T} = mod(atan(z, hypot(p.x, p.y) - r_torus), T(2π)) in θ


_in_torr_r_tube(p::CylindricalPoint, r_torus::Real, r_tube::Real) = hypot(p.r - r_torus, p.z) <= r_tube
_in_torr_r_tube(p::CylindricalPoint, r_torus::Real, r_tube::AbstractInterval) = hypot(p.r - r_torus, p.z) in r_tube

_in_torr_θ(p::CylindricalPoint{T}, r_torus::Real, θ::AbstractInterval) where {T} = mod(atan(z, p.r - r_torus), T(2π)) in θ


in(p::AbstractCoordinatePoint, t::Torus{<:Any, <:Any, <:Any, Nothing, Nothing}) =
    _in_torr_r_tube(p, t.r_torus, t.r_tube)

in(p::AbstractCoordinatePoint, t::Torus{<:Any, <:Any, <:Any, <:AbstractInterval, Nothing}) =
    _in_torr_r_tube(p, t.r_torus, t.r_tube) && _in_φ(p, t.φ)

in(p::AbstractCoordinatePoint, t::Torus{<:Any, <:Any, <:Any, <:AbstractInterval, <:AbstractInterval}) =
    _in_torr_r_tube(p, t.r_torus, t.r_tube) && _in_φ(p, t.φ) && _in_torr_θ(p, t.r_torus, t.θ)


# plotting
