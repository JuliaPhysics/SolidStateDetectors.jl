function get_plot_points(rp::RegularPrismMantle{N,T}; kwargs...) where {N, T <: AbstractFloat}

    plot_points = Vector{CartesianPoint{T}}[]
    zMin::T = _left_linear_interval(rp.z)
    zMax::T = _right_linear_interval(rp.z)

    #horizontal polygons
    for z in [zMin, zMax]
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rp.r * cos(φ), rp.r * sin(φ), z) for φ in range(0,2π,length=N+1)]))
    end

    #vertical lines
    for φ in range(2π/N,2π,length=N)
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rp.r * cos(φ), rp.r * sin(φ), zMin), CartesianPoint{T}(rp.r * cos(φ), rp.r * sin(φ), zMax)]))
    end

    plot_points
end

function mesh(rp::RegularPrismMantle{N,T}; n = 30) where {N, T <: AbstractFloat}

    zMin::T, zMax::T = get_z_limits(rp)
    φ = range(0, 2π, length = N+1)
    z = range(zMin, zMax, length = 2)

    X::Array{T,2} = [rp.r*cos(φ_i) for φ_i in φ, j in 1:length(z)]
    Y::Array{T,2} = [rp.r*sin(φ_i) for φ_i in φ, j in 1:length(z)]
    Z::Array{T,2} = [j for i in 1:length(φ), j in z]

    Mesh(X, Y, Z)
end
