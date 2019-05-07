# @recipe function f(Vol::SSD.ConeMantle{T}, n_aux_lines =0) where T
#
#     rStart = Vol.cone.r_interval.left
#     rStop = Vol.cone.r_interval.right
#     φStart = Vol.cone.φ_interval.left
#     φStop = Vol.cone.φ_interval.right
#     zStart = Vol.cone.z_interval.left
#     zStop = Vol.cone.z_interval.right
#
#     orientation = Vol.cone.orientation
#
#     if orientation == :bottom_left || orientation == :top_right
#         @series begin
#             line_3d(rStop,rStart,φStart,φStart,zStart,zStop)
#         end
#         @series begin
#             line_3d(rStop,rStart,φStop,φStop,zStart,zStop)
#         end
#         for ia in 0:n_aux_lines
#             @series begin
#                 line_3d(rStop,rStart,φStart+ia*(φStop-φStart)/(n_aux_lines+1),φStart+ia*(φStop-φStart)/(n_aux_lines+1),zStart,zStop)
#             end
#         end
#     end
# end
