function get_plot_points(rp::RegularPolygon{N,T}; kwargs...) where {N, T <: AbstractFloat}

    plot_points = Vector{CartesianPoint{T}}[]

    rMin::T = _left_radial_interval(rp.r)
    rMax::T = _right_radial_interval(rp.r)

    for r in (rMin == 0 ? [rMax] : [rMin, rMax])
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * cos(φ), r * sin(φ), rp.z) for φ in range(0,2π,length=N+1)]))
    end

    plot_points
end

function mesh(rp::RegularPolygon{N,T}; n = 30) where {N, T <: AbstractFloat}

    rMin::T, rMax::T = get_r_limits(rp)
    φrange = range(0, 2π, length = N+1)
    sφrange = sin.(φrange)
    cφrange = cos.(φrange)
    r = range(rMin, rMax, length = 2)
    z = fill(rp.z, length(r))

    X::Array{T,2} = [r_j*cφ for cφ in cφrange, r_j in r]
    Y::Array{T,2} = [r_j*sφ for sφ in sφrange, r_j in r]
    Z::Array{T,2} = [z_j for i in 1:N+1, z_j in z]

    Mesh(X, Y, Z)
end
