
function paint!(pointtypes, potential, face::AbstractPlanarSurfacePrimitive, geometry, pot_value, grid::CartesianGrid)
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
function paint!(pointtypes, potential, face::AbstractPlanarSurfacePrimitive, geometry, pot_value, grid::CylindricalGrid) where {T}
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
    if proj == Val{:rz}()
        p_temp = face.points[findmax(broadcast(p -> abs(p[1]), face.points))[2]]
        φ1 = CylindricalPoint(p_temp)[2]
        φ2 = φ1 < π ? φ1+π : φ1-π
        # Here we have to might have to do some modiciation of φ1 and φ2
        # due to the possible periodicity of the grid.
        for i2 in t_idx_r2
            for i1 in t_idx_r1
                for φ in (φ1, φ2)
                    pt = CylindricalPoint{T}(t1[i1], φ, t2[i2])
                    if pt in geometry
                        i3 = searchsortednearest(ticks[2], φ) 
                        pointtypes[i1, i3, i2] = zero(PointType)
                        potential[i1, i3, i2] = pot_value
                    end
                end
            end
        end
    else
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
                    end
                end
            end
        end 
    end
    nothing
end
