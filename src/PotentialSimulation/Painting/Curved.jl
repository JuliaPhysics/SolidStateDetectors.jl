function paint!(pointtypes, potential, face::AbstractSurfacePrimitive{T}, geometry, pot_value, grid::CartesianGrid) where {T}
    t_idx_r1, t_idx_r2, t_idx_r3 = get_sub_ind_ranges(face, grid)
    ticks = (grid.axes[1].ticks, grid.axes[2].ticks, grid.axes[3].ticks)
    eX = CartesianVector{T}(1,0,0);
    eY = CartesianVector{T}(0,1,0);
    eZ = CartesianVector{T}(0,0,1);
    widths_ax1 = diff(get_extended_ticks(grid[1]))
    widths_ax2 = diff(get_extended_ticks(grid[2]))
    widths_ax3 = diff(get_extended_ticks(grid[3]))
    for i2 in t_idx_r2
        for i1 in t_idx_r1
            l = ConstructiveSolidGeometry.Line(CartesianPoint{T}(ticks[1][i1], ticks[2][i2], zero(T)), eZ)
            pts = ConstructiveSolidGeometry.intersection(face, l)
            for pt in pts
                i3 = searchsortednearest(ticks[3], pt[3])
                csgtol = abs(ticks[3][i3] - pt[3])
                Δw_max = max(
                    widths_ax1[i1], widths_ax1[i1+1],
                    widths_ax2[i2], widths_ax2[i2+1],
                    # widths_ax3[i3], widths_ax3[i3+1]
                )/2
                if csgtol > Δw_max
                    # csgtol = Δw_max
                end
                if in(pt, geometry, csgtol = csgtol)
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
                i2 = searchsortednearest(ticks[2], pt[2])
                csgtol = abs(ticks[2][i2] - pt[2])
                Δw_max = max(
                    widths_ax1[i1], widths_ax1[i1+1],
                    # widths_ax2[i2], widths_ax2[i2+1],
                    widths_ax3[i3], widths_ax3[i3+1]
                )/2
                if csgtol > Δw_max
                    # csgtol = Δw_max
                end
                if in(pt, geometry, csgtol = csgtol)
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
                i1 = searchsortednearest(ticks[1], pt[1])
                csgtol = abs(ticks[1][i1] - pt[1])
                Δw_max = max(
                    # widths_ax1[i1], widths_ax1[i1+1],
                    widths_ax2[i2], widths_ax2[i2+1],
                    widths_ax3[i3], widths_ax3[i3+1]
                )/2
                if csgtol > Δw_max
                    # csgtol = Δw_max
                end
                if in(pt, geometry, csgtol = csgtol)
                    pointtypes[i1, i2, i3] = zero(PointType)
                    potential[i1, i2, i3] = pot_value
                end
            end
        end
    end 
    nothing
end


function paint!(pointtypes, potential, face::AbstractSurfacePrimitive{T}, geometry, pot_value, grid::CylindricalGrid) where {T}
    t_idx_r1, t_idx_r2, t_idx_r3 = get_sub_ind_ranges(face, grid)
    ticks = (grid.axes[1].ticks, grid.axes[2].ticks, grid.axes[3].ticks)
    eZ = CartesianVector{T}(0,0,1);

    widths_ax1 = diff(get_extended_ticks(grid[1]))
    widths_ax2 = diff(get_extended_ticks(grid[2]))
    widths_ax3 = diff(get_extended_ticks(grid[3]))

    for i2 in eachindex(ticks[2])
        for i1 in eachindex(ticks[1])
            l = ConstructiveSolidGeometry.Line(CartesianPoint(CylindricalPoint{T}(ticks[1][i1], ticks[2][i2], zero(T))), eZ)
            pts = ConstructiveSolidGeometry.intersection(face, l)
            for pt in pts
                i3 = searchsortednearest(ticks[3], pt[3])
                csgtol = abs(ticks[3][i3] - pt[3])
                Δw_max = max(
                    widths_ax1[i1], widths_ax1[i1+1],
                    # widths_ax2[i2], widths_ax2[i2+1],
                )/4
                if csgtol > Δw_max
                    csgtol = Δw_max
                end
                if in(pt, geometry, csgtol = csgtol)
                    pointtypes[i1, i2, i3] = zero(PointType)
                    potential[i1, i2, i3] = pot_value
                end
            end
        end
    end 
    # For this we need a the function `intersection_with_φ_axis`...
    # This will be fun... I skip this for now. 

    # for i3 in t_idx_r3 # z
    #     for i1 in t_idx_r1 # r
    #         l = ConstructiveSolidGeometry.Line(CartesianPoint{T}(ticks[1][i1], zero(T), ticks[3][i3]), eY)
    #         pt = ConstructiveSolidGeometry.intersection(plane, l)
#             i2 = searchsortednearest(ticks[2], pt[2])
    #         if in(pt, geometry, csgtol = abs(ticks[2][i2] - pt[2]))
    #             pointtypes[i1, i2, i3] = zero(PointType)
    #             potential[i1, i2, i3] = pot_value
    #         end
    #     end
    # end 

    for i3 in t_idx_r3 # z;   Maybe switch loops so that the direction of `l` has to be calculated less times..
        o = CartesianPoint{T}(zero(T), zero(T), ticks[3][i3])
        for i2 in eachindex(ticks[2]) # φ
            dir = CartesianVector(CartesianPoint(CylindricalPoint{T}(one(T), ticks[2][i2], zero(T)))) 
            l = ConstructiveSolidGeometry.Line(o, dir) # dir should be normalized
            pts_car = ConstructiveSolidGeometry.intersection(face, l)
            for pt_car in pts_car
                pt_cyl = CylindricalPoint(pt_car)
                # Intersection must be in positive r direction of dir: pt_cyl[2] == ticks[2][i2]
                # Use `abs(pt_cyl[2] - ticks[2][i2]) < 0.1` to avoid rounding issues
                # if it differs, it would always differ by π = 3.141... -> 0.1 is fine
                i1 = searchsortednearest(ticks[1], pt_cyl[1])
                csgtol = abs(ticks[1][i1] - pt_cyl[1])
                Δw_max = max(
                    widths_ax3[i3], widths_ax3[i3+1],
                    # widths_ax2[i2], widths_ax2[i2+1],
                )/4
                if csgtol > Δw_max
                    csgtol = Δw_max
                end
                if abs(pt_cyl[2] - ticks[2][i2]) < 0.1 && in(pt_car, geometry, csgtol = csgtol)
                    pointtypes[i1, i2, i3] = zero(PointType)
                    potential[i1, i2, i3] = pot_value
                end
            end
        end
    end 
    nothing
end
