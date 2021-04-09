struct ConeMantle{T,TR,TP,TZ} <: AbstractSurfacePrimitive{T}
    r::TR #if not a Tuple, then ConeMantle is a Tube
    φ::TP
    z::TZ
    #if r is a Tuple, the first entry refers to the r at the bottom, the second one to the r at the top
    function ConeMantle( ::Type{T},
                   r::Union{T, Tuple{T,T}},
                   φ::Union{Nothing, <:AbstractInterval{T}},
                   z::Union{T, <:AbstractInterval{T}}) where {T, I<:AbstractInterval{T}}
        new{T,typeof(r),typeof(φ),typeof(z)}(r, φ, z)
    end
end

function ConeMantle(c::Cone{T}; rbot = 1, rtop = 1) where {T}
    r = rbot == rtop ? T(rbot) : (T(rbot), T(rtop))
    ConeMantle( T, r, c.φ, c.z)
end

function ConeMantle(t::Torus{T}; θ = π/2) where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    θ = T(mod(θ,2π))
    sθ, cθ = sincos(θ)
    if θ > T(0) && θ < T(π)
        rbot = t.r_torus + r_tubeMin*cθ
        rtop = t.r_torus + r_tubeMax*cθ
        zMin = r_tubeMin*sθ
        zMax = r_tubeMax*sθ
    elseif θ > T(π) && θ < T(2π)
        rtop = t.r_torus + r_tubeMin*cθ
        rbot = t.r_torus + r_tubeMax*cθ
        zMax = r_tubeMin*sθ
        zMin = r_tubeMax*sθ
    else
        @error "Cone Mantle not defined for torroidal cordinate θ = 0 or θ = π. Use Annulus"
    end
    r = rbot == rtop ? T(rbot) : (T(rbot), T(rtop))
    z = T(zMin)..T(zMax)
    ConeMantle( T, r, t.φ, z)
end

function ConeMantle(;rbot = 1, rtop = 0, φMin = 0, φMax = 2π, zMin = -1/2, zMax = 1/2)
    T = float(promote_type(typeof.((rbot, rtop, φMin, φMax, zMin, zMax))...))
    r = rbot == rtop ? T(rbot) : (T(rbot), T(rtop))
    φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
    z = zMax == -zMin ? T(zMax) : T(zMin)..T(zMax)
    ConeMantle( T, r, φ, z)
end
ConeMantle(rbot, rtop, φMin, φMax, zMin, zMax) = ConeMantle(;rbot = rbot, rtop = rtop, φMin = φMin, φMax = φMax, zMin = zMin, zMax = zMax)

function ConeMantle(rbot::R1, rtop::R2, height::H) where {R1<:Real, R2<:Real, H<:Real}
    T = float(promote_type(R1, R2, H))
    ConeMantle( T, (T(rbot), T(rtop)), nothing, T(height)/2)
end

get_r_at_z(c::ConeMantle{T}, z::Real) where {T} = _get_r_at_z(get_r_limits(c)..., c.z, z)

get_r_limits(c::ConeMantle{T, T, <:Any, <:Any}) where {T} = (T(c.r), T(c.r))
get_r_limits(c::ConeMantle{T, <:Tuple, <:Any, <:Any}) where {T} = c.r

get_φ_limits(c::ConeMantle{T, <:Any, Nothing, <:Any}) where {T} = (T(0), T(2π), true)
get_φ_limits(c::ConeMantle{T, <:Any, <:AbstractInterval, <:Any}) where {T} = (c.φ.left, c.φ.right, false)

get_z_limits(c::ConeMantle{T}) where {T} = (_left_linear_interval(c.z), _right_linear_interval(c.z))

in(p::AbstractCoordinatePoint, c::ConeMantle{<:Any, <:Any, Nothing, <:Any}) =
    _in_z(p, c.z) && _eq_cyl_r(p, get_r_at_z(c, p.z))

in(p::AbstractCoordinatePoint, c::ConeMantle{<:Any, <:Any, <:AbstractInterval, <:Any}) =
    _in_z(p, c.z) && _in_φ(p, c.φ) && _eq_cyl_r(p, get_r_at_z(c, p.z))

function sample(c::ConeMantle{T}, step::Real) where {T}
    φMin::T, φMax::T, _ = get_φ_limits(c)
    zMin::T, zMax::T = get_z_limits(c)
    samples = [
        CylindricalPoint{T}(get_r_at_z(c, z),φ,z)
        for z in zMin:step:zMax
        for φ in (get_r_at_z(c, z) == 0 ? φMin : φMin:step/get_r_at_z(c, z):φMax)
    ]
end

function sample(c::ConeMantle{T}, Nsamps::NTuple{3,Int}) where {T}
    φMin::T, φMax::T, _ = get_φ_limits(c)
    zMin::T, zMax::T = get_z_limits(c)
    samples = [
        CylindricalPoint{T}(get_r_at_z(c, z),φ,z)
        for z in (Nsamps[3] ≤ 1 ? zMin : range(zMin, zMax, length = Nsamps[3]))
        for φ in (Nsamps[2] ≤ 1 ? φMin : range(φMin, φMax, length = Nsamps[2]))
    ]
end

function sample(c::ConeMantle{T}, g::CylindricalTuple{T}) where {T}
    samples = [
        CylindricalPoint{T}(get_r_at_z(c, z),φ,z)
        for z in get_z_ticks(c, g)
        for φ in get_φ_ticks(c, g)
    ]
end

function sample(a::ConeMantle{T}, g::CartesianTuple{T}) where {T}
    R::T = _right_radial_interval(a.r)
    x_from_y::Vector{T} = sqrt.(R^2 .- filter(y -> abs(y) <= R, g.y).^2)
    samples = [
        CartesianPoint{T}(x,y,z)
        for x in _get_ticks(sort!(vcat(g.x, x_from_y, -x_from_y)), -R, R)
        for y in [-sqrt(R^2-x^2),sqrt(R^2-x^2)]    
        for z in _get_ticks(g.z, _left_linear_interval(a.z), _right_linear_interval(a.z))
        if isnothing(a.φ) || mod(atan(y, x), T(2π)) in a.φ
    ]
end
