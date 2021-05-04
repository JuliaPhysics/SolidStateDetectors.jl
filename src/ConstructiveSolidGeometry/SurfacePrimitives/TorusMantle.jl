struct TorusMantle{T,TP,TT} <: AbstractSurfacePrimitive{T}
    r_torus::T
    r_tube::T
    φ::TP
    θ::TT
    z::T
    function TorusMantle( ::Type{T},
                   r_torus::T,
                   r_tube::T,
                   φ::Union{Nothing, <:AbstractInterval{T}},
                   θ::Union{Nothing, <:AbstractInterval{T}},
                   z::T) where {T}
        new{T,typeof(φ),typeof(θ)}(r_torus, r_tube, φ, θ, z)
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

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, Nothing, Nothing}) =
    _isapprox_torr_r_tube(p, t.r_torus, t.r_tube, t.z)

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, <:AbstractInterval, Nothing}) =
    _isapprox_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_φ(p, t.φ)

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, Nothing, <:AbstractInterval}) =
    _isapprox_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_torr_θ(p, t.r_torus, t.θ, t.z)

in(p::AbstractCoordinatePoint, t::TorusMantle{<:Any, <:AbstractInterval, <:AbstractInterval}) =
    _isapprox_torr_r_tube(p, t.r_torus, t.r_tube, t.z) && _in_φ(p, t.φ) && _in_torr_θ(p, t.r_torus, t.θ, t.z)


get_φ_limits(t::TorusMantle{T, Nothing}) where {T} = (T(0), T(2π), true)
get_φ_limits(t::TorusMantle{T, <:AbstractInterval}) where {T} = (t.φ.left, t.φ.right, false)

get_θ_limits(t::TorusMantle{T, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
get_θ_limits(t::TorusMantle{T, <:Any, <:AbstractInterval}) where {T} = (t.θ.left, t.θ.right, false)


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


function _get_z_ticks(t::TorusMantle{T}, g::CylindricalTicksTuple{T}) where {T}
    z_from_r::Vector{T} = sqrt.(t.r_tube.^2 .- (filter(r -> abs(r - t.r_torus) < t.r_tube, g.r).- t.r_torus).^2)
    filter!(z -> t.r_tube^2 - (z - t.z)^2 >= 0,_get_ticks(sort!(vcat(g.z, t.z .- z_from_r, t.z .+ z_from_r)), t.z - t.r_tube, t.z + t.r_tube))
end

function sample(t::TorusMantle{T}, g::CylindricalTicksTuple{T})::Vector{CylindricalPoint{T}} where {T}
    samples::Vector{CylindricalPoint{T}} = [
            CylindricalPoint{T}(r,φ,z)
            for z in _get_z_ticks(t, g)
            for φ in get_φ_ticks(t, g)
            for r in (t.r_torus - sqrt(t.r_tube^2 - (z - t.z)^2), t.r_torus + sqrt(t.r_tube^2 - (z - t.z)^2))
            if t.θ === nothing || _in_angular_interval_closed(mod(atan(z - t.z, r - t.r_torus), T(2π)), t.θ)
        ]
end


function _get_z_ticks(t::TorusMantle{T}, g::CartesianTicksTuple{T}) where {T}
    z_from_x::Vector{T} = t.z .+ sqrt.(t.r_tube.^2 .- (filter(x -> abs(x - t.r_torus) < t.r_tube, g.x).- t.r_torus).^2)
    z_from_y::Vector{T} = t.z .+ sqrt.(t.r_tube.^2 .- (filter(y -> abs(y - t.r_torus) < t.r_tube, g.y).- t.r_torus).^2)
    filter!(z -> t.r_tube^2 - (z - t.z)^2 >= 0,_get_ticks(unique!(sort!(vcat(g.z, t.z .- z_from_x, t.z .+ z_from_x, t.z .- z_from_y, t.z .+ z_from_y))), t.z - t.r_tube, t.z + t.r_tube))
end

function _get_x_at_z(t::TorusMantle{T}, g::CartesianTicksTuple{T}, z::T) where {T}
    R::T = sqrt(t.r_tube^2 - (z - t.z)^2)
    xMax_from_y::Vector{T} = sqrt.((t.r_torus + R)^2 .- filter(y -> abs(y) <= t.r_torus + R, g.y).^2)
    xMin_from_y::Vector{T} = sqrt.((t.r_torus - R)^2 .- filter(y -> abs(y) <= t.r_torus - R, g.y).^2)
    _get_ticks(sort!(vcat(-xMax_from_y, -xMin_from_y, g.x, xMin_from_y, xMax_from_y)), -t.r_torus - t.r_tube, t.r_torus + t.r_tube)
end

function _get_y_at_z(t::TorusMantle{T}, x::T, z::T) where {T}
    R::T = sqrt(t.r_tube^2 - (z - t.z)^2)
    tmp::T = (t.r_torus + R)^2 - x^2
    tmp2::T = (t.r_torus - R)^2 - x^2
    if tmp < 0 
        ()
    elseif tmp2 < 0 
        (-sqrt(tmp), sqrt(tmp))
    else 
        (-sqrt(tmp), -sqrt(tmp2), sqrt(tmp2), sqrt(tmp))
    end
end

function sample(t::TorusMantle{T}, g::CartesianTicksTuple{T})::Vector{CartesianPoint{T}} where {T}
    samples::Vector{CartesianPoint{T}} = [
            CartesianPoint{T}(x,y,z)
            for z in _get_z_ticks(t, g)
            for x in _get_x_at_z(t, g, z)
            for y in _get_y_at_z(t, x, z)
            if (t.φ === nothing || mod(atan(y, x), T(2π)) in t.φ) && 
               (t.θ === nothing || _in_angular_interval_closed(mod(atan(z - t.z, hypot(x, y) - t.r_torus), T(2π)), t.θ))
        ]
end

