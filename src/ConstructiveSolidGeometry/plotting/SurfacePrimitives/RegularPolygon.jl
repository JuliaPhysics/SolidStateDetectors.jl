function get_plot_points(rp::RegularPolygon{N,T}; kwargs...) where {N, T <: AbstractFloat}

    plot_points = Vector{CartesianPoint{T}}[]

    rMin::T, rMax::T = get_r_limits(rp)
    φrange = range(0, 2π, length = N+1)
    scφrange = sincos.(φrange)

    for r in (rMin == 0 ? rMax : (rMin, rMax))
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * cφ, r * sφ, rp.z) for (sφ,cφ) in scφrange]))
    end

    plot_points
end

function mesh(rp::RegularPolygon{N,T}; n = 30) where {N, T <: AbstractFloat}

    φrange = range(0, 2π, length = N+1)
    scφrange = sincos.(φrange)
    r = get_r_limits(rp)

    X::Array{T,2} = [r_j*cφ for (_,cφ) in scφrange, r_j in r]
    Y::Array{T,2} = [r_j*sφ for (sφ,_) in scφrange, r_j in r]
    Z::Array{T,2} = [rp.z for i in 1:N+1, _ in r]

    Mesh(X, Y, Z)
end
