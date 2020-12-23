function get_plot_points(h::HexagonalPrism{T}) where {T <: AbstractFloat}
    
    plot_points = Vector{CartesianPoint{T}}[]
    
    rMin::T = _left_radial_interval(h.r)
    rMax::T = _right_radial_interval(h.r)
    zMin::T = _left_linear_interval(h.z)
    zMax::T = _right_linear_interval(h.z)
   
    for r in (rMin == 0 ? [rMax] : [rMin, rMax])
        
        #horizontal hexagons
        for z in [zMin, zMax]
            push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * sin(φ), r * cos(φ), z) for φ in 0:π/3:2π]))
        end
            
        #vertical lines
        for φ in 0:π/3:5π/3
            push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * sin(φ), r * cos(φ), zMin), CartesianPoint{T}(r * sin(φ), r * cos(φ), zMax)]))
        end
    end

    plot_points
end



# Old
# @recipe function f(hp::HexagonalPrism{T}) where {T <: SSDFloat}
#     label --> "HexagonalPrism"
#     @series begin
#         pts_top_outer = []
#         pts_bottom_outer = []
#         pts_top_inner = []
#         pts_bottom_inner = []

#         #find all vertices, this loop has been tested and works
#         for φ in deg2rad.(30:60:330) .- hp.rotZ
#             pt_top_outer = CartesianPoint{T}(hp.translate.x + hp.rOuter * cos(φ), hp.translate.y + hp.rOuter * sin(φ), hp.translate.z + hp.h/2)
#             push!(pts_top_outer, pt_top_outer)
#             pt_top_inner = CartesianPoint{T}(hp.translate.x + hp.rInner * cos(φ), hp.translate.y + hp.rInner * sin(φ), hp.translate.z + hp.h/2)
#             push!(pts_top_inner, pt_top_inner)
#             pt_bottom_outer = CartesianPoint{T}(hp.translate.x + hp.rOuter * cos(φ), hp.translate.y + hp.rOuter * sin(φ), hp.translate.z -hp.h/2)
#             push!(pts_bottom_outer, pt_bottom_outer)
#             pt_bottom_inner = CartesianPoint{T}(hp.translate.x + hp.rInner * cos(φ), hp.translate.y + hp.rInner * sin(φ), hp.translate.z -hp.h/2)
#             push!(pts_bottom_inner, pt_bottom_inner)
#         end

#         #create Linesegments connecting the vertices
#         lines = LineSegment{T, 3, :cartesian}[]
#         N = length(pts_top_outer)
#         for i in 1:N
#             push!(lines, LineSegment(pts_top_outer[i%N+1], pts_top_outer[i])) #top outer hexagon
#             push!(lines, LineSegment(pts_bottom_outer[i%N+1], pts_bottom_outer[i])) #bottom outer hexagon
#             push!(lines, LineSegment(pts_bottom_outer[i], pts_top_outer[i])) #lines connecting outer hexagons
#             if hp.rInner > 0
#                 push!(lines, LineSegment(pts_bottom_inner[i%N+1], pts_bottom_inner[i])) #bottom inner hexagon
#                 push!(lines, LineSegment(pts_top_inner[i%N+1], pts_top_inner[i])) #top inner hexagon
#                 push!(lines, LineSegment(pts_bottom_inner[i], pts_top_inner[i])) #lines connecting inner hexagons
#                 push!(lines, LineSegment(pts_top_outer[i], pts_top_inner[i])) #lines connecting hexagons
#                 push!(lines, LineSegment(pts_bottom_outer[i], pts_bottom_inner[i])) #lines connecting hexagons
#             end
#         end
#         lines
#     end
# end