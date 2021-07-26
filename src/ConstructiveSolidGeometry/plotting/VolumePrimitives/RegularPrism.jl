function get_plot_points(rp::RegularPrism{N,T}; kwargs...) where {N, T <: AbstractFloat}
    plot_points = Vector{CartesianPoint{T}}[]
    zMin::T, zMax::T = _linear_endpoints(rp.hZ)

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
