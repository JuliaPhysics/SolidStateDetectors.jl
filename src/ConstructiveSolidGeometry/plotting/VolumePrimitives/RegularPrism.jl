function get_plot_points(rp::RegularPrism{N,T}; kwargs...) where {N, T <: AbstractFloat}
    
    plot_points = Vector{CartesianPoint{T}}[]
    
    rMin::T = _left_radial_interval(rp.r)
    rMax::T = _right_radial_interval(rp.r)
    zMin::T = _left_linear_interval(rp.z)
    zMax::T = _right_linear_interval(rp.z)
   
    for r in (rMin == 0 ? [rMax] : [rMin, rMax])
        
        #horizontal polygons
        for z in [zMin, zMax]
            push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * cos(φ), r * sin(φ), z) for φ in range(0,2π,length=N+1)]))
        end
            
        #vertical lines
        for φ in range(2π/N,2π,length=N)
            push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * cos(φ), r * sin(φ), zMin), CartesianPoint{T}(r * cos(φ), r * sin(φ), zMax)]))
        end
    end

    plot_points
end