struct CylindricalAnnulus{T,TR,TP,TZ} <: AbstractSurfacePrimitive{T}
    r::TR
    φ::TP
    z::TZ
    function CylindricalAnnulus( ::Type{T},
                   r::Union{T, <:AbstractInterval{T}},
                   φ::Union{Nothing, <:AbstractInterval{T}},
                   z::T) where {T}
        new{T,typeof(r),typeof(φ),T}(r, φ, z)
    end
end

#Constructors
CylindricalAnnulus(c::Cone{T}; z = 0) where {T} = CylindricalAnnulus(T, get_r_at_z(c,z), c.φ, T(z))

function CylindricalAnnulus(t::Torus{T}; θ = 0) where {T}
    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    if θ == T(0)
        rMin = t.r_torus + r_tubeMin
        rMax = t.r_torus + r_tubeMax
    elseif θ == T(π)
        rMin = t.r_torus - r_tubeMax
        rMax = t.r_torus - r_tubeMin
    else
        @error "CylindricalAnnulus not defined for torroidal cordinate θ ≠ 0 and θ ≠ π. Use ConeMantle"
    end
    r = rMin == 0 ? T(rMax) : T(rMin)..T(rMax)
    CylindricalAnnulus( T, r, t.φ, T(0))
end


function CylindricalAnnulus(; rMin = 0, rMax = 1, φMin = 0, φMax = 2π, z = 0)
    T = float(promote_type(typeof.((rMin, rMax, φMin, φMax, z))...))
    r = rMin == 0 ? T(rMax) : T(rMin)..T(rMax)
    φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
    CylindricalAnnulus(T, r, φ, T(z))
end

CylindricalAnnulus(rMin, rMax, φMin, φMax, z) = CylindricalAnnulus(;rMin = rMin, rMax = rMax, φMin = φMin, φMax = φMax, z = z)

function CylindricalAnnulus(r::Real, z::Real)
    T = float(promote_type(typeof.((r, z))...))
    CylindricalAnnulus(T, T(r), nothing, T(z), nothing)
end

get_r_limits(a::CylindricalAnnulus{T, <:Union{T, AbstractInterval{T}}, <:Any}) where {T} =
    (_left_radial_interval(a.r), _right_radial_interval(a.r))

get_φ_limits(a::CylindricalAnnulus{T, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
get_φ_limits(a::CylindricalAnnulus{T, <:Any, <:AbstractInterval}) where {T} = (a.φ.left, a.φ.right, false)

in(p::AbstractCoordinatePoint, a::CylindricalAnnulus{T, <:Any, Nothing}) where {T} = _eq_z(p, a.z) && _in_cyl_r(p, a.r)

in(p::AbstractCoordinatePoint, a::CylindricalAnnulus{T, <:Any, <:AbstractInterval}) where {T} = _eq_z(p, a.z) && _in_φ(p, a.φ) && _in_cyl_r(p, a.r)

function distance_to_surface(point::AbstractCoordinatePoint, a::CylindricalAnnulus{T, <:Any, Nothing})::T where {T}
    point = CylindricalPoint(point)
    rMin::T, rMax::T = get_r_limits(a)
    _in_cyl_r(point, a.r) ? abs(point.z - a.z) : hypot(point.z - a.z, min(abs(point.r - rMin), abs(point.r - rMax)))
end

function distance_to_surface(point::AbstractCoordinatePoint, a::CylindricalAnnulus{T, <:Any, <:AbstractInterval})::T where {T}
    point = CylindricalPoint(point)
    rMin::T, rMax::T = get_r_limits(a)
    φMin::T, φMax::T, _ = get_φ_limits(a)
    Δz = abs(point.z - a.z)
    if _in_φ(point, a.φ)
        d = _in_cyl_r(point, a.r) ? Δz : hypot(Δz, min(abs(point.r - rMin), abs(point.r - rMax)))
    else
        ΔφMin = mod(point.φ - φMin, T(2π))
        ΔφMax = mod(point.φ - φMax, T(2π))
        Δφ = min(min(ΔφMin, T(2π) - ΔφMin), min(ΔφMax, T(2π) - ΔφMax))
        y, x = point.r .* sincos(Δφ)
        d = if x < rMin
            sqrt((rMin - x)^2 + y^2 +  Δz^2)
        elseif x > rMax
            sqrt((rMax - x)^2 + y^2 +  Δz^2)
        else
            hypot(y, Δz)
        end
    end
    d
end

function sample(a::CylindricalAnnulus{T}, step::Quantity{<:Real, Unitful.𝐋}) where {T}
    samples = CylindricalPoint{T}[]
    rMin::T, rMax::T = get_r_limits(a)
    φMin::T, φMax::T, _ = get_φ_limits(a)
    step = T(ustrip(uconvert(u"m", step)))
    for r in rMin:step:rMax
        if r == 0
            push!(samples, CylindricalPoint{T}(0,0,a.z))
        else
            for φ in φMin:step/r:φMax
                push!(samples, CylindricalPoint{T}(r,φ,a.z))
            end
        end
    end
    samples
end

#plotting
function get_plot_points(a::CylindricalAnnulus{T}; n = 30) where {T <: AbstractFloat}

    plot_points = Vector{CartesianPoint{T}}[]

    rMin::T, rMax::T = get_r_limits(a)
    φMin::T, φMax::T, φ_is_full_2π::Bool = get_φ_limits(a)
    φrange = range(φMin, φMax, length = n + 1)

    #circle(s)
    for r in [rMin, rMax]
        if r == 0 continue end
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * cos(φ), r * sin(φ), a.z) for φ in φrange]))
    end

    #for incomplete φ: lines of cross-sections
    if !φ_is_full_2π
        for φ in [φMin, φMax]
            push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rMin * cos(φ), rMin * sin(φ), a.z), CartesianPoint{T}(rMax * cos(φ), rMax * sin(φ), a.z)]))
        end
    end
    plot_points
end

function mesh(a::CylindricalAnnulus{T}; n = 30) where {T <: AbstractFloat}

    rMin::T, rMax::T = get_r_limits(a)
    φMin::T, φMax::T, _ = get_φ_limits(a)

    φrange = range(φMin, φMax, length = n+1)
    sφrange = sin.(φrange)
    cφrange = cos.(φrange)
    r = range(rMin, rMax, length = 2)
    z = fill(a.z, length(r))

    X::Array{T,2} = [r_j*cφ for cφ in cφrange, r_j in r]
    Y::Array{T,2} = [r_j*sφ for sφ in sφrange, r_j in r]
    Z::Array{T,2} = [z_j for i in 1:n+1, z_j in z]

    Mesh(X, Y, Z)
end
