
function paint!(pointtypes, potential, face::AbstractCurvedSurfacePrimitive{T}, geometry, pot_value, grid::CartesianGrid) where {T}
    t_idx_r1, t_idx_r2, t_idx_r3 = get_sub_ind_ranges(face, grid)
    ticks = (grid.axes[1].ticks, grid.axes[2].ticks, grid.axes[3].ticks)
    eX = CartesianVector{T}(1,0,0);
    eY = CartesianVector{T}(0,1,0);
    eZ = CartesianVector{T}(0,0,1);
    for i2 in t_idx_r2
        for i1 in t_idx_r1
            l = ConstructiveSolidGeometry.Line(CartesianPoint{T}(ticks[1][i1], ticks[2][i2], zero(T)), eZ)
            pts = ConstructiveSolidGeometry.intersection(face, l)
            for pt in pts
                if pt in geometry
                    i3 = searchsortednearest(ticks[3], pt[3])
                    pointtypes[i1, i2, i3] = zero(PointType)
                    potential[i1, i2, i3] = pot_value
                end
            end
        end
    end 
    for i3 in t_idx_r3
        for i1 in t_idx_r1
            l = ConstructiveSolidGeometry.Line(CartesianPoint{T}(ticks[1][i1], zero(T), ticks[3][i3]), eY)
            pts = ConstructiveSolidGeometry.intersection(face, l)
            for pt in pts
                if pt in geometry
                    i2 = searchsortednearest(ticks[2], pt[2])
                    pointtypes[i1, i2, i3] = zero(PointType)
                    potential[i1, i2, i3] = pot_value
                end
            end
        end
    end 
    for i3 in t_idx_r3
        for i2 in t_idx_r2
            l = ConstructiveSolidGeometry.Line(CartesianPoint{T}(zero(T), ticks[2][i2], ticks[3][i3]), eX)
            pts = ConstructiveSolidGeometry.intersection(face, l)
            for pt in pts
                if pt in geometry
                    i1 = searchsortednearest(ticks[1], pt[1])
                    pointtypes[i1, i2, i3] = zero(PointType)
                    potential[i1, i2, i3] = pot_value
                end
            end
        end
    end 
    nothing
end


# function paint!(pointtypes, potential, face::ConstructiveSolidGeometry.CylinderMantle, geometry, pot_value, grid::CartesianGrid)
#     ticks = TicksTuple(grid)
#     t_idx_r1, t_idx_r2, proj = ConstructiveSolidGeometry.get_2d_grid_ticks_and_proj(face, ticks)
#     t1, t2 = if proj == Val{:xy}() 
#         ticks[1], ticks[2]
#     elseif proj == Val{:xz}() 
#         ticks[1], ticks[3]
#     else
#         ticks[2], ticks[3]
#     end
#     for i2 in t_idx_r2
#         for i1 in t_idx_r1
#             pt1, pt2 = ConstructiveSolidGeometry.evaluate(face, t1[i1], t2[i2], proj)
#             for pt in (pt1, pt2)
#                 if pt in geometry
#                     if proj == Val{:xy}() 
#                         i3 = searchsortednearest(ticks[3], pt[3])
#                         pointtypes[i1, i2, i3] = zero(PointType)
#                         potential[i1, i2, i3] = pot_value
#                     elseif proj == Val{:xz}() 
#                         i3 = searchsortednearest(ticks[2], pt[2])
#                         pointtypes[i1, i3, i2] = zero(PointType)
#                         potential[i1, i3, i2] = pot_value
#                     else
#                         i3 = searchsortednearest(ticks[1], pt[1])
#                         pointtypes[i3, i1, i2] = zero(PointType)
#                         potential[i3, i1, i2] = pot_value
#                     end
#                 end
#             end
#         end
#     end 
#     nothing
# end