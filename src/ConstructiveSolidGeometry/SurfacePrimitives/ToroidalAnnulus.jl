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

#plotting
function get_plot_points(t::ToroidalAnnulus{T}; n = 30) where {T <: AbstractFloat}

    plot_points = Vector{CartesianPoint{T}}[]

    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    θMin::T, θMax::T, θ_is_full_2π::Bool = get_θ_limits(t)
    θrange = range(θMin, θMax, length = n + 1)

    #circle(s)
    for r_tube in [r_tubeMin, r_tubeMax]
        if r_tube == 0 continue end
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}((t.r_torus+r_tube*cos(θ))cos(t.φ), (t.r_torus+r_tube*cos(θ))sin(t.φ), r_tube*sin(θ)) for θ in θrange]))
    end

    #for incomplete θ: lines of cross-sections
    if !θ_is_full_2π
        for θ in [θMin, θMax]
            push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}((t.r_torus+r_tubeMin*cos(θ))cos(t.φ), (t.r_torus+r_tubeMin*cos(θ))sin(t.φ), r_tubeMin*sin(θ)), CartesianPoint{T}((t.r_torus+r_tubeMax*cos(θ))cos(t.φ), (t.r_torus+r_tubeMax*cos(θ))sin(t.φ), r_tubeMax*sin(θ))]))
        end
    end
    plot_points
end

function mesh(t::ToroidalAnnulus{T}; n = 30) where {T <: AbstractFloat}

    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    r_tube = range(r_tubeMin, r_tubeMax, length = 2)
    sφ, cφ = sincos(t.φ)
    θrange = range(θMin, θMax, length = n + 1)
    sθrange = sin.(θrange)
    cθrange = cos.(θrange)

    X::Array{T,2} = [(t.r_torus+r*cθ)*cφ for cθ in cθrange, r in r_tube]
    Y::Array{T,2} = [(t.r_torus+r*cθ)*sφ for cθ in cθrange, r in r_tube]
    Z::Array{T,2} = [r*sθ for sθ in sθrange, r in r_tube]

    Mesh(X, Y, Z)
end
