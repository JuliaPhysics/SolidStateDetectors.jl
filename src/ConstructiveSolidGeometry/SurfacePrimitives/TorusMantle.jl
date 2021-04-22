struct TorusMantle{T,TR,TB,TP,TT,TZ} <: AbstractSurfacePrimitive{T}
    r_torus::TR
    r_tube::TB
    φ::TP
    θ::TT
    z::TZ
    function TorusMantle( ::Type{T},
                   r_torus::T,
                   r_tube::T,
                   φ::Union{Nothing, <:AbstractInterval{T}},
                   θ::Union{Nothing, <:AbstractInterval{T}},
                   z::T) where {T}
        new{T,T,T,typeof(φ),typeof(θ),T}(r_torus, r_tube, φ, θ, z)
    end
end

#Constructors
TorusMantle(t::Torus{T}; r_tube = 1) where {T} = TorusMantle( T, t.r_torus, T(r_tube), t.φ, t.θ, t.z)

function TorusMantle(;r_torus = 1, r_tube = 1, φMin = 0, φMax = 2π, θMin = 0, θMax = 2π, z = 0)
    T = float(promote_type(typeof.((r_torus, r_tube, φMin, φMax, θMin, θMax, z))...))
    φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
    θ = mod(T(θMax) - T(θMin), T(2π)) == 0 ? nothing : T(θMin)..T(θMax)
    TorusMantle( T, T(r_torus), T(r_tube), φ, θ, T(z))
end
TorusMantle(r_torus, r_tube, φMin, φMax, θMin, θMax, z) = TorusMantle(;r_torus = r_torus, r_tube = r_tube, φMin = φMin, φMax = φMax, θMin = θMin, θMax = θMax, z = z)

function get_surface_vector(t::TorusMantle{T}, point::AbstractCoordinatePoint)::CartesianVector{T} where {T}
    pcy = CylindricalPoint(point)
    r = pcy.r - t.r_torus
    sφ::T, cφ::T = sincos(pcy.φ)
    CartesianVector{T}(r*cφ, r*sφ, pcy.z - t.z)
end

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, <:Any, <:Any, Nothing, Nothing}) =
    _eq_torr_r_tube(p, t.r_torus, t.r_tube, t.z)

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, <:Any, <:Any, <:AbstractInterval, Nothing}) =
    _eq_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_φ(p, t.φ)

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, <:Any, <:Any, Nothing, <:AbstractInterval}) =
    _eq_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_torr_θ(p, t.r_torus, t.θ, t.z)

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, <:Any, <:Any, <:AbstractInterval, <:AbstractInterval}) =
    _eq_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_φ(p, t.φ) && _in_torr_θ(p, t.r_torus, t.θ, t.z)


get_φ_limits(t::TorusMantle{T, <:Any, <:Any, Nothing, <:Any}) where {T} = (T(0), T(2π), true)
get_φ_limits(t::TorusMantle{T, <:Any, <:Any, <:AbstractInterval, <:Any}) where {T} = (t.φ.left, t.φ.right, false)

get_θ_limits(t::TorusMantle{T, <:Any, <:Any, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
get_θ_limits(t::TorusMantle{T, <:Any, <:Any, <:Any, <:AbstractInterval}) where {T} = (t.θ.left, t.θ.right, false)


function sample(t::TorusMantle{T}, step::Real) where {T}
    φMin::T, φMax::T, _ = get_φ_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    samples = [
        CylindricalPoint{T}(get_r_at_θ(t,θ),φ,t.r_tube*sin(θ)+t.z)
        for φ in (t.r_tube == 0 ? φMin : φMin:step/t.r_tube:φMax)
        for θ in (t.r_tube == 0 ? θMin : θMin:step/t.r_tube:θMax)
    ]
end

function sample(t::TorusMantle{T}, Nsamps::NTuple{3,Int}) where {T}
    φMin::T, φMax::T, _ = get_φ_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    samples = [
        CylindricalPoint{T}(get_r_at_θ(t,θ),φ,t.r_tube*sin(θ)+t.z)
        for φ in (Nsamps[2] ≤ 1 ? φMin : range(φMin, φMax, length = Nsamps[2]))
        for θ in (Nsamps[3] ≤ 1 ? θMin : range(θMin, θMax, length = Nsamps[3]))
    ]
end
