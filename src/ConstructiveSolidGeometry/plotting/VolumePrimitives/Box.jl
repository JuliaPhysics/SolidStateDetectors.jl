function get_plot_points(b::Box{T}) where {T <: AbstractFloat}
    
    plot_points = Vector{CartesianPoint{T}}[]
    
    xMin::T = _left_linear_interval(b.x)
    xMax::T = _right_linear_interval(b.x)
    yMin::T = _left_linear_interval(b.y)
    yMax::T = _right_linear_interval(b.y)
    zMin::T = _left_linear_interval(b.z)
    zMax::T = _right_linear_interval(b.z)
    
    #horizontal squares
    for z in [zMin, zMax]
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(xMin, yMin, z), CartesianPoint{T}(xMin, yMax, z),
                            CartesianPoint{T}(xMax, yMax, z), CartesianPoint{T}(xMax, yMin, z), 
                            CartesianPoint{T}(xMin, yMin, z)]))
    end
    
    #vertical lines
    for x in [xMin, xMax]
        for y in [yMin, yMax]
            push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(x,y,zMin), CartesianPoint{T}(x,y,zMax)]))
        end
    end
        
    plot_points
end


# Old CSG
# @recipe function f(cb::Box{T}) where {T <: SSDFloat}
#     label-->"Box"
#     ls = LineSegments(cb)
#     @series begin
#         ls
#     end
# end