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
# plotting

function get_plot_points(t::TorusMantle{T}; n = 30) where {T <: AbstractFloat}
    plot_points = Vector{CartesianPoint{T}}[]
    φMin::T, φMax::T, φ_is_full_2π::Bool = get_φ_limits(t)
    θMin::T, θMax::T, θ_is_full_2π::Bool = get_θ_limits(t)
    sθMin, cθMin = sincos(θMin)

    r1 = T(t.r_torus + t.r_tube*cθMin)
    z1 = T(t.r_tube*sθMin)
    θ2 = θ_is_full_2π ? π : θMax
    sθ2, cθ2 = sincos(θ2)
    r2 = T(t.r_torus + t.r_tube*cθ2)
    z2 = T(t.r_tube*sθ2)

    append!(plot_points, get_plot_points(CylindricalAnnulus(T,r1..r1,t.φ,z1), n = n))
    append!(plot_points, get_plot_points(CylindricalAnnulus(T,r2..r2,t.φ,z2), n = n))
    for φ in (φ_is_full_2π ? [φMin] : [φMin, φMax])
        append!(plot_points, get_plot_points(ToroidalAnnulus(T,t.r_torus,t.r_tube..t.r_tube,φ,t.θ), n = n))
    end
    plot_points
end

function mesh(t::TorusMantle{T}; n = 30) where {T <: AbstractFloat}
    φMin::T, φMax::T, _ = get_φ_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    θrange = range(θMin, θMax, length = n + 1)
    sθrange = sin.(θrange)
    cθrange = cos.(θrange)
    φrange = range(φMin, φMax, length = n + 1)
    sφrange = sin.(φrange)
    cφrange = cos.(φrange)

    X = [(t.r_torus + t.r_tube*cθ)*cφ for cφ in cφrange, cθ in cθrange]
    Y = [(t.r_torus + t.r_tube*cθ)*sφ for sφ in sφrange, cθ in cθrange]
    Z = [t.r_tube*sθ for i in 1:n+1, sθ in sθrange]
    Mesh(X, Y, Z)
end
