
function paint!(pointtypes, potential, face::ConstructiveSolidGeometry.Polygon, geometry, pot_value, grid::CartesianGrid)
    ticks = TicksTuple(grid)
    t_idx_r1, t_idx_r2, proj = ConstructiveSolidGeometry.get_2d_grid_ticks_and_proj(face, ticks)
    t1, t2 = if proj == Val{:xy}() 
        ticks[1], ticks[2]
    elseif proj == Val{:xz}() 
        ticks[1], ticks[3]
    else
        ticks[2], ticks[3]
    end
    plane = ConstructiveSolidGeometry.Plane(face)
    for i2 in t_idx_r2
        for i1 in t_idx_r1
            pt = ConstructiveSolidGeometry.evaluate(plane, t1[i1], t2[i2], proj)
            if pt in geometry
                if proj == Val{:xy}() 
                    i3 = searchsortednearest(ticks[3], pt[3])
                    pointtypes[i1, i2, i3] = zero(PointType)
                    potential[i1, i2, i3] = pot_value
                elseif proj == Val{:xz}() 
                    i3 = searchsortednearest(ticks[2], pt[2])
                    pointtypes[i1, i3, i2] = zero(PointType)
                    potential[i1, i3, i2] = pot_value
                else
                    i3 = searchsortednearest(ticks[1], pt[1])
                    pointtypes[i3, i1, i2] = zero(PointType)
                    potential[i3, i1, i2] = pot_value
                end
            end
        end
    end 
    nothing
end
function paint!(pointtypes, potential, face::ConstructiveSolidGeometry.Polygon, geometry, pot_value, grid::CylindricalGrid)
    ticks = TicksTuple(grid)
    t_idx_r1, t_idx_r2, proj = ConstructiveSolidGeometry.get_2d_grid_ticks_and_proj(face, ticks)
    t1, t2 = if proj == Val{:rφ}() 
        ticks[1], ticks[2]
    elseif proj == Val{:φz}() 
        ticks[2], ticks[3]
    else
        ticks[1], ticks[3]
    end
    plane = ConstructiveSolidGeometry.Plane(face)
    for i2 in t_idx_r2
        for i1 in t_idx_r1
            pt = ConstructiveSolidGeometry.evaluate(plane, t1[i1], t2[i2], proj)
            if pt in geometry
                if proj == Val{:rφ}() 
                    i3 = searchsortednearest(ticks[3], pt[3])
                    pointtypes[i1, i2, i3] = zero(PointType)
                    potential[i1, i2, i3] = pot_value
                elseif proj == Val{:φz}() 
                    i3 = searchsortednearest(ticks[1], pt[1])
                    pointtypes[i3, i1, i2] = zero(PointType)
                    potential[i3, i1, i2] = pot_value
                else # proj == Val{:rz}() # Do we need this case? I will leave it for now. 
                    i3 = searchsortednearest(ticks[2], pt[2])
                    pointtypes[i1, i3, i2] = zero(PointType)
                    potential[i1, i3, i2] = pot_value
                end
            end
        end
    end 
    nothing
end
