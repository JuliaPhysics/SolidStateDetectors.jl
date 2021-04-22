function get_plot_points(t::TorusMantle{T}; n = 30) where {T <: AbstractFloat}
    plot_points = Vector{CartesianPoint{T}}[]
    φMin::T, φMax::T, φ_is_full_2π::Bool = get_φ_limits(t)
    θMin::T, θMax::T, θ_is_full_2π::Bool = get_θ_limits(t)
    sθMin, cθMin = sincos(θMin)

    r1 = T(t.r_torus + t.r_tube*cθMin)
    z1 = T(t.r_tube*sθMin) + t.z
    θ2 = θ_is_full_2π ? π : θMax
    sθ2, cθ2 = sincos(θ2)
    r2 = T(t.r_torus + t.r_tube*cθ2)
    z2 = T(t.r_tube*sθ2) + t.z

    append!(plot_points, get_plot_points(CylindricalAnnulus(T,r1..r1,t.φ,z1), n = n))
    append!(plot_points, get_plot_points(CylindricalAnnulus(T,r2..r2,t.φ,z2), n = n))
    for φ in (φ_is_full_2π ? [φMin] : [φMin, φMax])
        append!(plot_points, get_plot_points(ToroidalAnnulus(T,t.r_torus,t.r_tube..t.r_tube,φ,t.θ,t.z), n = n))
    end
    plot_points
end

function mesh(t::TorusMantle{T}; n = 30) where {T <: AbstractFloat}
    φMin::T, φMax::T, _ = get_φ_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)

    fφ = (φMax - φMin)/(2π)
    nφ = Int(ceil(n*fφ))

    fθ = (θMax - θMin)/(2π)
    nθ = Int(ceil(n*fθ))

    θrange = range(θMin, θMax, length = nθ + 1)
    sθrange = sin.(θrange)
    cθrange = cos.(θrange)
    φrange = range(φMin, φMax, length = nφ + 1)
    sφrange = sin.(φrange)
    cφrange = cos.(φrange)

    X = [(t.r_torus + t.r_tube*cθ)*cφ for cφ in cφrange, cθ in cθrange]
    Y = [(t.r_torus + t.r_tube*cθ)*sφ for sφ in sφrange, cθ in cθrange]
    Z = [t.r_tube*sθ + t.z for i in 1:nφ+1, sθ in sθrange]
    Mesh(X, Y, Z)
end
