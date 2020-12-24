function get_plot_points(s::Sphere{T}) where {T <: AbstractFloat}
    
    plot_points = Vector{CartesianPoint{T}}[]
    
    rMin::T = _left_radial_interval(s.r)
    rMax::T = _right_radial_interval(s.r)
    φrange = range(0, 2π, length = 36)
    
    for r in (rMin == 0 ? [rMax] : [rMin, rMax])
        for φ in range(0, 2π, length = 11)
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * sin(θ) * cos(φ), r * sin(θ) * sin(φ), r * cos(θ)) for θ in φrange]))
        end
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * cos(φ), r * sin(φ), 0) for φ in φrange]))
    end

    plot_points
end

# Old
# @recipe function f(s::Sphere{T};) where {T <: SSDFloat}
#     @series begin
#         pts = []
#         label --> ""
#         for φ in range(T(0), length = 36, stop = T(2π))
#             pt = s.org + CartesianVector{T}( s.r * cos(φ), s.r * sin(φ), 0  )
#             push!(pts, pt)
#         end
#         lines = LineSegment{T, 3, :cartesian}[]
#         for i in 1:length(pts)-1
#             push!(lines, LineSegment(pts[i+1], pts[i]))
#         end
#         lines
#     end
#     @series begin
#         label --> ""
#         pts = []
#         for φ in range(T(0), length = 36, stop = T(2π))
#             pt = s.org + CartesianVector{T}( s.r * cos(φ), 0, s.r * sin(φ) )
#             push!(pts, pt)
#         end
#         lines = LineSegment{T, 3, :cartesian}[]
#         for i in 1:length(pts)-1
#             push!(lines, LineSegment(pts[i+1], pts[i]))
#         end
#         lines
#     end
#     @series begin
#         pts = []
#         label --> ""      
#         for φ in range(T(0), length = 36, stop = T(2π))
#             pt = s.org + CartesianVector{T}( 0, s.r * cos(φ), s.r * sin(φ) )
#             push!(pts, pt)
#         end
#         lines = LineSegment{T, 3, :cartesian}[]
#         for i in 1:length(pts)-1
#             push!(lines, LineSegment(pts[i+1], pts[i]))
#         end
#         lines
#     end
# end