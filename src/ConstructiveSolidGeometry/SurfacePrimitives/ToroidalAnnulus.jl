struct ToroidalAnnulus{T,TR,TB,TP,TT} <: AbstractSurfacePrimitive{T}
    r_torus::TR
    r_tube::TB
    φ::TP
    θ::TT
    function ToroidalAnnulus( ::Type{T},
                   r_torus::T,
                   r_tube::Union{T, <:AbstractInterval{T}},
                   φ::T,
                   θ::Union{Nothing, <:AbstractInterval{T}}) where {T}
        new{T,T,typeof(r_tube),T,typeof(θ)}(r_torus, r_tube, φ, θ)
    end
end

#Constructors
ToroidalAnnulus(t::Torus{T}; φ = 0) where {T} = ToroidalAnnulus( T, t.r_torus, t.r_tube, T(mod(φ,2π)), t.θ)

function ToroidalAnnulus(;r_torus = 1, r_tubeMin = 0, r_tubeMax = 1, φ = 0, θMin = 0, θMax = 2π)
    T = float(promote_type(typeof.((r_torus, r_tubeMin, r_tubeMax, φ, θMin, θMax))...))
    r_tube = r_tubeMin == 0 ? T(r_tubeMax) : T(r_tubeMin)..T(r_tubeMax)
    θ = mod(T(θMax) - T(θMin), T(2π)) == 0 ? nothing : T(θMin)..T(θMax)
    ToroidalAnnulus( T, T(r_torus), r_tube, T(φ), θ)
end
ToroidalAnnulus(r_torus, r_tubeMin, r_tubeMax, φ, θMin, θMax) = ToroidalAnnulus(;r_torus = r_torus, r_tubeMin = r_tubeMin, r_tubeMax = r_tubeMax, φ = φ, θMin = θMin, θMax = θMax)

in(p::AbstractCoordinatePoint, t::ToroidalAnnulus{<:Any, <:Any, <:Any, <:Any, Nothing}) =
    _in_torr_r_tube(p, t.r_torus, t.r_tube) && _eq_φ(p, t.φ)

in(p::AbstractCoordinatePoint, t::ToroidalAnnulus{<:Any, <:Any, <:Any, <:Any, <:AbstractInterval}) =
    _in_torr_r_tube(p, t.r_torus, t.r_tube) && _eq_φ(p, t.φ) && _in_torr_θ(p, t.r_torus, t.θ)

get_r_tube_limits(t::ToroidalAnnulus{T}) where {T} =
    (_left_radial_interval(t.r_tube),_right_radial_interval(t.r_tube))

get_θ_limits(t::ToroidalAnnulus{T, <:Any, <:Any, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
get_θ_limits(t::ToroidalAnnulus{T, <:Any, <:Any, <:Any, <:AbstractInterval}) where {T} = (t.θ.left, t.θ.right, false)

function sample(t::ToroidalAnnulus{T}, step::Real) where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    samples = [
        CylindricalPoint{T}(t.r_torus+r_tube*cos(θ),t.φ,r_tube*sin(θ))
        for r_tube in r_tubeMin:step:r_tubeMax
        for θ in (r_tube == 0 ? θMin : θMin:step/r_tube:θMax)
    ]
end

function sample(t::ToroidalAnnulus{T}, Nsamps::NTuple{3,Int}) where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    samples = [
        CylindricalPoint{T}(t.r_torus+r_tube*cos(θ),t.φ,r_tube*sin(θ))
        for r_tube in (Nsamps[1] ≤ 1 ? r_tubeMin : range(r_tubeMin, r_tubeMax, length = Nsamps[1]))
        for θ in (Nsamps[3] ≤ 1 ? θMin : range(θMin, θMax, length = Nsamps[3]))
    ]
end
