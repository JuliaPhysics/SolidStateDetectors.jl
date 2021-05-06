struct Arc{T,TP,TH} <: AbstractLinePrimitive{T}
    r::T
    center::TP
    α::TH
    function Arc( ::Type{T},
                   r::T,
                   center::PlanarPoint{T},
                   α::Union{Nothing, <:AbstractInterval{T}}) where {T}
        new{T,typeof(center),typeof(α)}(r, center, α)
    end
end

function Arc(; r = 0, center = PlanarPoint(0,0), αMin = 0, αMax = 2π)
    T = float(promote_type(typeof.((r, αMin, αMax))..., eltype(center)))
    α = mod(T(αMax) - T(αMin), T(2π)) == 0 ? nothing : T(αMin)..T(αMax)
    Arc(T, T(r), PlanarPoint{T}(center), α)
end

Arc(r, center, αMin, αMax) = Arc(; r = r, center = center, αMin = αMin, αMax = αMax)

Circle(r::T, center::PlanarPoint{T}) where {T} = Arc(T, r, center, nothing)
Circle(a::Arc{T}) where {T} = Arc(T, a.r, a.center, nothing)

get_α_at_u_v(a::Arc{T}, u::Real, v::Real) where {T} = mod(atan(v - a.center.v, u - a.center.u), 2π) #u,v are planar coordinates

get_α_limits(a::Arc{T, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
get_α_limits(a::Arc{T, <:Any, <:AbstractInterval}) where {T} = (a.α.left, a.α.right, false)


function sample(a::Arc{T}, step::AbstractFloat) where {T}
    αMin::T, αMax::T, _ = get_α_limits(a)
    samples = [PlanarPoint{T}(a.r*cos(α)+a.center.u,a.r*sin(α)+a.center.v) for α in (a.r == 0 ? αMin : αMin:step/a.r:αMax)]
end

function sample(a::Arc{T}, Nsamps::Int) where {T}
    αMin::T, αMax::T, _ = get_α_limits(a)
    samples = [PlanarPoint{T}(a.r*cos(α)+a.center.u,a.r*sin(α)+a.center.v) for α in range(αMin, αMax, length = Nsamps)]
end

function distance_to_line(point::PlanarPoint{T}, a::Arc{T, <:Any, Nothing}) where {T}
    p_t = point - a.center
    r = hypot(p_t.u, p_t.v)
    abs(a.r - r)
end

function distance_to_line(point::PlanarPoint{T}, a::Arc{T, <:Any, <:AbstractInterval}) where {T}
    αMin::T, αMax::T, _ = get_α_limits(a)
    p_t = point - a.center
    r = hypot(p_t.u, p_t.v)
    α = atan(p_t.v, p_t.u)
    if _in_angular_interval_closed(α, a.α)
        return abs(a.r - r)
    else
        αNear = Δ_φ(T(α),αMin) ≤ Δ_φ(T(α),αMax) ? αMin : αMax
        sαNear, cαNear = sincos(αNear)
        return norm(p_t - PlanarPoint{T}(a.r*cαNear, a.r*sαNear))
    end
end
