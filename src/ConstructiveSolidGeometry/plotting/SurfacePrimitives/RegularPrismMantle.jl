function get_plot_points(rp::RegularPrismMantle{N,T}; kwargs...) where {N, T <: AbstractFloat}

    plot_points = Vector{CartesianPoint{T}}[]
    zMin::T, zMax::T = get_z_limits(rp)

    φrange = range(0, 2π, length = N+1)
    scφrange = sincos.(φrange)

    #horizontal polygons
    for z in (zMin, zMax)
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rp.r * cφ, rp.r * sφ, z) for (sφ,cφ) in scφrange]))
    end

    #vertical lines
    for (sφ,cφ) in scφrange[1:end-1]
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rp.r * cφ, rp.r * sφ, zMin), CartesianPoint{T}(rp.r * cφ, rp.r * sφ, zMax)]))
    end

    plot_points
end

function mesh(rp::RegularPrismMantle{N,T}; n = 30) where {N, T <: AbstractFloat}

    φrange = range(0, 2π, length = N+1)
    scφrange = sincos.(φrange)
    z = get_z_limits(rp)

    X::Array{T,2} = [rp.r*cφ for (_,cφ) in scφrange, _ in z]
    Y::Array{T,2} = [rp.r*sφ for (sφ,_) in scφrange, _ in z]
    Z::Array{T,2} = [j for i in φrange, j in z]

    Mesh(X, Y, Z)
end
