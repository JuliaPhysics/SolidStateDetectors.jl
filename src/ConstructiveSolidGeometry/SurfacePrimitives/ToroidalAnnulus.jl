struct ToroidalAnnulus{T,TB,TT} <: AbstractSurfacePrimitive{T}
    r_torus::T
    r_tube::TB
    φ::T
    θ::TT
    z::T
    function ToroidalAnnulus( ::Type{T},
                   r_torus::T,
                   r_tube::Union{T, <:AbstractInterval{T}},
                   φ::T,
                   θ::Union{Nothing, <:AbstractInterval{T}},
                   z::T) where {T}
        new{T,typeof(r_tube),typeof(θ)}(r_torus, r_tube, φ, θ, z)
    end
end

#Constructors
ToroidalAnnulus(t::Torus{T}; φ = 0) where {T} = ToroidalAnnulus( T, t.r_torus, t.r_tube, T(mod(φ,2π)), t.θ, t.z)

function ToroidalAnnulus(;r_torus = 1, r_tubeMin = 0, r_tubeMax = 1, φ = 0, θMin = 0, θMax = 2π, z = 0)
    T = float(promote_type(typeof.((r_torus, r_tubeMin, r_tubeMax, φ, θMin, θMax, z))...))
    r_tube = r_tubeMin == 0 ? T(r_tubeMax) : T(r_tubeMin)..T(r_tubeMax)
    θ = mod(T(θMax) - T(θMin), T(2π)) == 0 ? nothing : T(θMin)..T(θMax)
    ToroidalAnnulus( T, T(r_torus), r_tube, T(φ), θ, T(z))
end

ToroidalAnnulus(r_torus, r_tubeMin, r_tubeMax, φ, θMin, θMax, z) = ToroidalAnnulus(;r_torus = r_torus, r_tubeMin = r_tubeMin, r_tubeMax = r_tubeMax, φ = φ, θMin = θMin, θMax = θMax, z = z)

in(p::AbstractCoordinatePoint, t::ToroidalAnnulus{<:Any, <:Any, Nothing}) =
    _in_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _isapprox_φ(p, t.φ)

in(p::AbstractCoordinatePoint, t::ToroidalAnnulus{<:Any, <:Any, <:AbstractInterval}) =
    _in_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _isapprox_φ(p, t.φ) && _in_torr_θ(p, t.r_torus, t.θ, t.z)

get_r_tube_limits(t::ToroidalAnnulus{T}) where {T} =
    (_left_radial_interval(t.r_tube),_right_radial_interval(t.r_tube))

get_θ_limits(t::ToroidalAnnulus{T, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
get_θ_limits(t::ToroidalAnnulus{T, <:Any, <:AbstractInterval}) where {T} = (t.θ.left, t.θ.right, false)
get_z_limits(t::ToroidalAnnulus{T}) where {T} = (t.z - _right_radial_interval(t.r_tube), t.z + _right_radial_interval(t.r_tube))

function sample(t::ToroidalAnnulus{T}, step::Real) where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    samples = [
        CylindricalPoint{T}(t.r_torus+r_tube*cos(θ),t.φ,r_tube*sin(θ)+t.z)
        for r_tube in r_tubeMin:step:r_tubeMax
        for θ in (r_tube == 0 ? θMin : θMin:step/r_tube:θMax)
    ]
end

function sample(t::ToroidalAnnulus{T}, Nsamps::NTuple{3,Int}) where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    samples = [
        CylindricalPoint{T}(t.r_torus+r_tube*cos(θ),t.φ,r_tube*sin(θ)+t.z)
        for r_tube in (Nsamps[1] ≤ 1 ? r_tubeMin : range(r_tubeMin, r_tubeMax, length = Nsamps[1]))
        for θ in (Nsamps[3] ≤ 1 ? θMin : range(θMin, θMax, length = Nsamps[3]))
    ]
end

function _get_z_at_r(t::ToroidalAnnulus{T,T}, g::CylindricalTicksTuple{T}, r::T) where {T}
    tmp::T = t.r_tube^2 - (r - t.r_torus)^2 
    if tmp < 0 return (t.z,) end
    _get_ticks(g.z, t.z - sqrt(tmp), t.z + sqrt(tmp))
end

function _get_z_at_r(t::ToroidalAnnulus{T,<:AbstractInterval{T}}, g::CylindricalTicksTuple{T}, r::T) where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    tmp::T = r_tubeMax^2 - (r - t.r_torus)^2 
    if tmp < 0 return (t.z,) end
    tmp2::T = r_tubeMin^2 - (r - t.r_torus)^2 
    if tmp2 < 0 return _get_ticks(g.z, t.z - sqrt(tmp), t.z + sqrt(tmp)) end
    vcat(_get_ticks(g.z, t.z - sqrt(tmp), t.z - sqrt(tmp2)),_get_ticks(g.z, t.z + sqrt(tmp2), t.z + sqrt(tmp)))
end
    

function sample(t::ToroidalAnnulus{T}, g::CylindricalTicksTuple{T})::Vector{CylindricalPoint{T}} where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    samples::Vector{CylindricalPoint{T}} = [
            CylindricalPoint{T}(r,t.φ,z)
            for φ in _get_ticks(g.φ, t.φ, t.φ) # only sample if t.φ is within the grid bounds
            for r in _get_ticks(g.r, t.r_torus - r_tubeMax, t.r_torus + r_tubeMax)
            for z in _get_z_at_r(t, g, r)
            if t.θ === nothing || _in_angular_interval_closed(mod(atan(z - t.z, r - t.r_torus), T(2π)), t.θ)
        ]
end

function sample(t::ToroidalAnnulus{T}, g::CartesianTicksTuple{T})::Vector{CartesianPoint{T}} where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    sφ::T, cφ::T = sincos(t.φ)
    r_ticks = unique!(sort!(vcat((cφ == 0 ? [] : g.x./cφ), (sφ == 0 ? [] : g.y./sφ))))
    samples::Vector{CartesianPoint{T}} = [
        CartesianPoint{T}(r*cφ,r*sφ,z)
        for z in get_z_ticks(t, g)
        for r in _get_ticks(r_ticks, t.r_torus - r_tubeMax, t.r_torus + r_tubeMax)
        if (r_tubeMin <= hypot(r - t.r_torus, z - t.z) <= r_tubeMax) &&
           (t.θ === nothing || _in_angular_interval_closed(mod(atan(z - t.z, r - t.r_torus), T(2π)), t.θ))
    ]
end