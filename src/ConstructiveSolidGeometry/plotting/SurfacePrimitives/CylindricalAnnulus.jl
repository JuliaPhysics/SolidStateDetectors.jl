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
    f = (φMax - φMin)/(2π)
    n = Int(ceil(n*f))
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
