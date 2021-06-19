
function paint!(pointtypes, potential, face::AbstractPlanarSurfacePrimitive{T}, geometry, pot_value, grid::CartesianGrid) where {T}
    t_idx_r1, t_idx_r2, t_idx_r3 = get_sub_ind_ranges(face, grid)
    ticks = (grid.axes[1].ticks, grid.axes[2].ticks, grid.axes[3].ticks)
    plane = ConstructiveSolidGeometry.Plane(face)
    eX = CartesianVector{T}(1,0,0);
    eY = CartesianVector{T}(0,1,0);
    eZ = CartesianVector{T}(0,0,1);
    for i2 in t_idx_r2
        for i1 in t_idx_r1
            l = ConstructiveSolidGeometry.Line(CartesianPoint{T}(ticks[1][i1], ticks[2][i2], zero(T)), eZ)
            pt = ConstructiveSolidGeometry.intersection(plane, l)
            if pt in geometry
                i3 = searchsortednearest(ticks[3], pt[3])
                pointtypes[i1, i2, i3] = zero(PointType)
                potential[i1, i2, i3] = pot_value
            end
        end
    end 
    for i3 in t_idx_r3
        for i1 in t_idx_r1
            l = ConstructiveSolidGeometry.Line(CartesianPoint{T}(ticks[1][i1], zero(T), ticks[3][i3]), eY)
            pt = ConstructiveSolidGeometry.intersection(plane, l)
            if pt in geometry
                i2 = searchsortednearest(ticks[2], pt[2])
                pointtypes[i1, i2, i3] = zero(PointType)
                potential[i1, i2, i3] = pot_value
            end
        end
    end 
    for i3 in t_idx_r3
        for i2 in t_idx_r2
            l = ConstructiveSolidGeometry.Line(CartesianPoint{T}(zero(T), ticks[2][i2], ticks[3][i3]), eX)
            pt = ConstructiveSolidGeometry.intersection(plane, l)
            if pt in geometry
                i1 = searchsortednearest(ticks[1], pt[1])
                pointtypes[i1, i2, i3] = zero(PointType)
                potential[i1, i2, i3] = pot_value
            end
        end
    end 
    nothing
end


function paint!(pointtypes, potential, face::AbstractPlanarSurfacePrimitive{T}, geometry, pot_value, grid::CylindricalGrid) where {T}
    t_idx_r1, t_idx_r2, t_idx_r3 = get_sub_ind_ranges(face, grid)
    ticks = (grid.axes[1].ticks, grid.axes[2].ticks, grid.axes[3].ticks)
    plane = ConstructiveSolidGeometry.Plane(face)
    eZ = CartesianVector{T}(0,0,1);
    for i2 in eachindex(ticks[2])
        for i1 in eachindex(ticks[1])
            l = ConstructiveSolidGeometry.Line(CartesianPoint{T}(ticks[1][i1], ticks[2][i2], zero(T)), eZ)
            pt = ConstructiveSolidGeometry.intersection(plane, l)
            if pt in geometry
                i3 = searchsortednearest(ticks[3], pt[3])
                pointtypes[i1, i2, i3] = zero(PointType)
                potential[i1, i2, i3] = pot_value
            end
        end
    end 
    # for i3 in t_idx_r3 # z
    #     for i1 in t_idx_r1 # r
        # For this we need a new method for `intersection(plane, Circle)`
    #         l = ConstructiveSolidGeometry.Line(CartesianPoint{T}(ticks[1][i1], zero(T), ticks[3][i3]), eY)
    #         pt = ConstructiveSolidGeometry.intersection(plane, l)
    #         if pt in geometry
    #             i2 = searchsortednearest(ticks[2], pt[2])
    #             pointtypes[i1, i2, i3] = zero(PointType)
    #             potential[i1, i2, i3] = pot_value
    #         end
    #     end
    # end 
    for i3 in t_idx_r3 # z;   Maybe switch loops so that the direction of `l` has to be calculated less times..
        o = CartesianPoint{T}(zero(T), zero(T), ticks[3][i3])
        for i2 in eachindex(ticks[2]) # φ
            l = ConstructiveSolidGeometry.Line(o, CartesianVector(CartesianPoint(CylindricalPoint{T}(one(T), ticks[2][i2], zero(T)))))
            pt = CylindricalPoint(ConstructiveSolidGeometry.intersection(plane, l))
            if pt in geometry
                i1 = searchsortednearest(ticks[1], pt[1])
                pointtypes[i1, i2, i3] = zero(PointType)
                potential[i1, i2, i3] = pot_value
            end
        end
    end 
    nothing
end
