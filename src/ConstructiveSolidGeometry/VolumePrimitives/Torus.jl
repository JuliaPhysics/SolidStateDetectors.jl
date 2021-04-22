struct Torus{T,TR,TB,TP,TT,TZ} <: AbstractVolumePrimitive{T}
    r_torus::TR
    r_tube::TB
    φ::TP
    θ::TT
    z::TZ
    function Torus( ::Type{T},
                   r_torus::T,
                   r_tube::Union{T, <:AbstractInterval{T}},
                   φ::Union{Nothing, <:AbstractInterval{T}},
                   θ::Union{Nothing, <:AbstractInterval{T}},
                   z::T) where {T}
        new{T,T,typeof(r_tube),typeof(φ),typeof(θ),T}(r_torus, r_tube, φ, θ, z)
    end
end

#Constructors
function Torus(;r_torus = 1, r_tubeMin = 0, r_tubeMax = 1, φMin = 0, φMax = 2π, θMin = 0, θMax = 2π, z = 0)
    T = float(promote_type(typeof.((r_torus, r_tubeMin, r_tubeMax, φMin, φMax, θMin, θMax, z))...))
    r_tube = r_tubeMin == 0 ? T(r_tubeMax) : T(r_tubeMin)..T(r_tubeMax)
    φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
    θ = mod(T(θMax) - T(θMin), T(2π)) == 0 ? nothing : T(θMin)..T(θMax)
    Torus( T, T(r_torus), r_tube, φ, θ, T(z))
end

Torus(r_torus, r_tubeMin, r_tubeMax, φMin, φMax, θMin, θMax, z) = Torus(;r_torus = r_torus, r_tubeMin = r_tubeMin, r_tubeMax = r_tubeMax, φMin = φMin, φMax = φMax, θMin = θMin, θMax = θMax, z = z)

function Torus(r_torus::R1, r_tube::R2, z::TZ) where {R1<:Real, R2<:Real, TZ<:Real}
    T = float(promote_type(R1, R2, TZ))
    Torus( T, T(r_torus), T(r_tube), nothing, nothing, T(z))
end

function RoundChamfer(r_torus::R1, r_tube::R2, z::TZ) where {R1<:Real, R2<:Real, TZ<:Real}
    T = float(promote_type(R1, R2, TZ))
    Torus( T, T(r_torus), T(r_tube), nothing, T(0)..T(π/2), T(z))
end

in(p::AbstractCoordinatePoint, t::Torus{<:Any, <:Any, <:Any, Nothing, Nothing}) =
    _in_torr_r_tube(p, t.r_torus, t.r_tube, t.z)

in(p::AbstractCoordinatePoint, t::Torus{<:Any, <:Any, <:Any, <:AbstractInterval, Nothing}) =
    _in_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_φ(p, t.φ)

in(p::AbstractCoordinatePoint, t::Torus{<:Any, <:Any, <:Any, Nothing, <:AbstractInterval}) =
    _in_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_torr_θ(p, t.r_torus, t.θ, t.z)

in(p::AbstractCoordinatePoint, t::Torus{<:Any, <:Any, <:Any, <:AbstractInterval, <:AbstractInterval}) =
    _in_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_φ(p, t.φ) && _in_torr_θ(p, t.r_torus, t.θ, t.z)

get_r_tube_limits(t::Torus{T}) where {T} = (_left_radial_interval(t.r_tube),_right_radial_interval(t.r_tube))

get_φ_limits(t::Torus{T, <:Any, <:Any, Nothing, <:Any}) where {T} = (T(0), T(2π), true)
get_φ_limits(t::Torus{T, <:Any, <:Any, <:AbstractInterval, <:Any}) where {T} = (t.φ.left, t.φ.right, false)

get_θ_limits(t::Torus{T, <:Any, <:Any, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
get_θ_limits(t::Torus{T, <:Any, <:Any, <:Any, <:AbstractInterval}) where {T} = (t.θ.left, t.θ.right, false)

function _get_decomposed_surfaces(t::Torus{T}) where {T}
    surfaces = AbstractSurfacePrimitive[]
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    for r_tube in [r_tubeMin, r_tubeMax]
        if r_tube == 0 continue end
        push!(surfaces, TorusMantle(t, r_tube = r_tube))
    end
    unique(surfaces)
end

get_decomposed_surfaces(t::Torus{T, T, <:Any, Nothing, Nothing}) where {T} = _get_decomposed_surfaces(t)

function get_decomposed_surfaces(t::Torus{T, T, <:Any, <:AbstractInterval, Nothing}) where {T}
    φMin::T, φMax::T, _ = get_φ_limits(t)
    surfaces = _get_decomposed_surfaces(t)
    for φ in [φMin, φMax]
        push!(surfaces, ToroidalAnnulus(t, φ = φ))
    end
    unique(surfaces)
end

function get_decomposed_surfaces(t::Torus{T, T, <:Any, Nothing, <:AbstractInterval}) where {T}
    θMin::T, θMax::T, _ = get_θ_limits(t)
    surfaces = _get_decomposed_surfaces(t)
    for θ in [θMin, θMax]
        θ in [T(0),T(π)] ? push!(surfaces, CylindricalAnnulus(t, θ = θ)) : push!(surfaces, ConeMantle(t, θ = θ))
    end
    unique(surfaces)
end

function get_decomposed_surfaces(t::Torus{T, T, <:Any, <:AbstractInterval, <:AbstractInterval}) where {T}
    θMin::T, θMax::T, _ = get_θ_limits(t)
    φMin::T, φMax::T, _ = get_φ_limits(t)
    surfaces = _get_decomposed_surfaces(t)
    for θ in [θMin, θMax]
        θ in [T(0),T(π)] ? push!(surfaces, CylindricalAnnulus(t, θ = θ)) : push!(surfaces, ConeMantle(t, θ = θ))
    end
    for φ in [φMin, φMax]
        push!(surfaces, ToroidalAnnulus(t, φ = φ))
    end
    unique(surfaces)
end

function sample(t::Torus{T}, step::Real) where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    φMin::T, φMax::T, _ = get_φ_limits(t)
    samples = [
        CylindricalPoint{T}(t.r_torus+r_tube*cos(θ),φ,r_tube*sin(θ))
        for r_tube in r_tubeMin:step:r_tubeMax
        for θ in (r_tube == 0 ? θMin : θMin:step/r_tube:θMax)
        for φ in (r_tube == 0 ? φMin : φMin:step/r_tube:φMax)
    ]
end

function sample(t::Torus{T}, Nsamps::NTuple{3,Int} = (2,5,3)) where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    φMin::T, φMax::T, _ = get_φ_limits(t)
    samples = [
        CylindricalPoint{T}(t.r_torus+r_tube*cos(θ),φ,r_tube*sin(θ))
        for r_tube in (Nsamps[1] ≤ 1 ? r_tubeMin : range(r_tubeMin, r_tubeMax, length = Nsamps[1]))
        for θ in (Nsamps[3] ≤ 1 ? θMin : range(θMin, θMax, length = Nsamps[3]))
        for φ in (Nsamps[2] ≤ 1 ? φMin : range(φMin, φMax, length = Nsamps[2]))
    ]
end
