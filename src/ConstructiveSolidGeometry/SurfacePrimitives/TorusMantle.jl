struct TorusMantle{T,TR,TB,TP,TT} <: AbstractSurfacePrimitive{T}
    r_torus::TR
    r_tube::TB
    φ::TP
    θ::TT
    function TorusMantle( ::Type{T},
                   r_torus::T,
                   r_tube::T,
                   φ::Union{Nothing, <:AbstractInterval{T}},
                   θ::Union{Nothing, <:AbstractInterval{T}}) where {T}
        new{T,T,T,typeof(φ),typeof(θ)}(r_torus, r_tube, φ, θ)
    end
end

#Constructors
TorusMantle(t::Torus{T}; r_tube = 1) where {T} = TorusMantle( T, t.r_torus, T(r_tube), t.φ, t.θ)

function TorusMantle(;r_torus = 1, r_tube = 1, φMin = 0, φMax = 2π, θMin = 0, θMax = 2π)
    T = float(promote_type(typeof.((r_torus, r_tube, φMin, φMax, θMin, θMax))...))
    φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
    θ = mod(T(θMax) - T(θMin), T(2π)) == 0 ? nothing : T(θMin)..T(θMax)
    TorusMantle( T, T(r_torus), T(r_tube), φ, θ)
end
TorusMantle(r_torus, r_tube, φMin, φMax, θMin, θMax) = TorusMantle(;r_torus = r_torus, r_tube = r_tube, φMin = φMin, φMax = φMax, θMin = θMin, θMax = θMax)

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, <:Any, <:Any, Nothing, Nothing}) =
    _eq_torr_r_tube(p, t.r_torus, t.r_tube)

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, <:Any, <:Any, <:AbstractInterval, Nothing}) =
    _eq_torr_r_tube(p, t.r_torus, t.r_tube) && _in_φ(p, t.φ)

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, <:Any, <:Any, Nothing, <:AbstractInterval}) =
    _eq_torr_r_tube(p, t.r_torus, t.r_tube) && _in_torr_θ(p, t.r_torus, t.θ)

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, <:Any, <:Any, <:AbstractInterval, <:AbstractInterval}) =
    _eq_torr_r_tube(p, t.r_torus, t.r_tube) && _in_φ(p, t.φ) && _in_torr_θ(p, t.r_torus, t.θ)


get_φ_limits(t::TorusMantle{T, <:Any, <:Any, Nothing, <:Any}) where {T} = (T(0), T(2π), true)
get_φ_limits(t::TorusMantle{T, <:Any, <:Any, <:AbstractInterval, <:Any}) where {T} = (t.φ.left, t.φ.right, false)

get_θ_limits(t::TorusMantle{T, <:Any, <:Any, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
get_θ_limits(t::TorusMantle{T, <:Any, <:Any, <:Any, <:AbstractInterval}) where {T} = (t.θ.left, t.θ.right, false)
