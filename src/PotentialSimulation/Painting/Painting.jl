function get_sub_ind_ranges(a::Union{
        <:Tuple{AbstractVector{<:CartesianPoint{T}}, CartesianTicksTuple{T}},
        <:Tuple{AbstractVector{<:CylindricalPoint{T}}, CylindricalTicksTuple{T}}
        }) where {T}
    # For Cylindrical Coordinates this only works for z right now
    # There, we will only use t_idx_range_ax3 in the painting (for now...)
    pts, t = a[1], a[2]
    ax1, ax2, ax3 = getindex.(pts, 1), getindex.(pts, 2), getindex.(pts, 3)
    # ± 2 because due to merging close ticks of an axis of a grid,
    # the grid points which should "belong" to corners of a primitive
    # might be a bit outside of the intervals defined by the points 
    # of the `extreme_points`-method.
    # ± 1 might also be enough -> testing...
    ax1_min_idx = searchsortedfirst(t[1], minimum(ax1)) - 2
    ax2_min_idx = searchsortedfirst(t[2], minimum(ax2)) - 2
    ax3_min_idx = searchsortedfirst(t[3], minimum(ax3)) - 2
    ax1_max_idx = searchsortedfirst(t[1], maximum(ax1)) + 2 
    ax2_max_idx = searchsortedfirst(t[2], maximum(ax2)) + 2 
    ax3_max_idx = searchsortedfirst(t[3], maximum(ax3)) + 2 
    ls = (length(t[1]), length(t[2]), length(t[3]))
    if ax1_min_idx < 1 ax1_min_idx = 1 end
    if ax2_min_idx < 1 ax2_min_idx = 1 end
    if ax3_min_idx < 1 ax3_min_idx = 1 end
    if ax1_max_idx > ls[1] ax1_max_idx = ls[1] end
    if ax2_max_idx > ls[2] ax2_max_idx = ls[2] end
    if ax3_max_idx > ls[3] ax3_max_idx = ls[3] end
    t_idx_range_ax1 = ax1_min_idx:ax1_max_idx
    t_idx_range_ax2 = ax2_min_idx:ax2_max_idx
    t_idx_range_ax3 = ax3_min_idx:ax3_max_idx
    t_idx_range_ax1, t_idx_range_ax2, t_idx_range_ax3
end

get_sub_ind_ranges(p::ConstructiveSolidGeometry.AbstractSurfacePrimitive{T}, grid::CartesianGrid3D{T}) where {N,T} = 
    get_sub_ind_ranges((ConstructiveSolidGeometry.extreme_points(p), TicksTuple(grid)))
get_sub_ind_ranges(p::ConstructiveSolidGeometry.AbstractSurfacePrimitive{T}, grid::CylindricalGrid{T}) where {N,T} = 
    get_sub_ind_ranges((CylindricalPoint.(ConstructiveSolidGeometry.extreme_points(p)), TicksTuple(grid)))

function paint!(point_types, potential, face::AbstractSurfacePrimitive{T}, geometry, pot_value, grid::CartesianGrid3D{T}) where {T}
    t_idx_r1, t_idx_r2, t_idx_r3 = get_sub_ind_ranges(face, grid)
    ticks = (grid.axes[1].ticks, grid.axes[2].ticks, grid.axes[3].ticks)
    eX = CartesianVector{T}(1,0,0);
    eY = CartesianVector{T}(0,1,0);
    eZ = CartesianVector{T}(0,0,1);
    widths_ax1 = diff(get_extended_ticks(grid[1]))
    widths_ax2 = diff(get_extended_ticks(grid[2]))
    widths_ax3 = diff(get_extended_ticks(grid[3]))
    Δw_max_factor = T(1e-2)
    #= 
        Δw_max_factor is chosen by trying out different values for it. 
        This value seems to be okay if the grid is not to unevenly spaced. 
        But this can be secured via the `max_ratio`- and the `the max_tick_distance`-keywords.
        The critical parameter is `csgtol` inside the `in`-method, which is currently defined through 
        the widths of the voxel of the grid point next to the calculated intersection point.
        It would be probably better to turn this into an `NTuple{3,T}` and pass the widths of the voxel instead
        of one single value (`csgtol`). However, inside the `in` methods this would have to be transformed 
        due to rotations and maybe the shape of the primitive. This might not be straight forward.
        But this can be investigated and certainly improved in the future.
    =#
    for i2 in t_idx_r2
        for i1 in t_idx_r1
            l = ConstructiveSolidGeometry.Line(CartesianPoint{T}(ticks[1][i1], ticks[2][i2], zero(T)), eZ)
            pts = ConstructiveSolidGeometry.intersection(face, l)
            for pt in pts
                i3 = searchsortednearest(ticks[3], pt[3])
                csgtol = abs(ticks[3][i3] - pt[3])
                Δw_max = Δw_max_factor * max(
                    widths_ax1[i1], widths_ax1[i1+1],
                    widths_ax2[i2], widths_ax2[i2+1],
                    # widths_ax3[i3], widths_ax3[i3+1]
                ) 
                if csgtol < Δw_max
                    csgtol = Δw_max
                end
                if in(pt, geometry, csgtol = csgtol)
                    point_types[i1, i2, i3] = zero(PointType)
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
                Δw_max = Δw_max_factor * max(
                    widths_ax1[i1], widths_ax1[i1+1],
                    # widths_ax2[i2], widths_ax2[i2+1],
                    widths_ax3[i3], widths_ax3[i3+1]
                )
                if csgtol < Δw_max
                    csgtol = Δw_max
                end
                if in(pt, geometry, csgtol = csgtol)
                    point_types[i1, i2, i3] = zero(PointType)
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
                Δw_max = Δw_max_factor * max(
                    # widths_ax1[i1], widths_ax1[i1+1],
                    widths_ax2[i2], widths_ax2[i2+1],
                    widths_ax3[i3], widths_ax3[i3+1]
                )
                if csgtol < Δw_max
                    csgtol = Δw_max
                end
                if in(pt, geometry, csgtol = csgtol)
                    point_types[i1, i2, i3] = zero(PointType)
                    potential[i1, i2, i3] = pot_value
                end
            end
        end
    end 
    nothing
end


function paint!(point_types, potential, face::AbstractSurfacePrimitive{T}, geometry, pot_value, grid::CylindricalGrid) where {T}
    t_idx_r1, t_idx_r2, t_idx_r3 = get_sub_ind_ranges(face, grid)
    ticks = (grid.axes[1].ticks, grid.axes[2].ticks, grid.axes[3].ticks)
    eZ = CartesianVector{T}(0,0,1);

    widths_ax1 = diff(get_extended_ticks(grid[1]))
    widths_ax2 = diff(get_extended_ticks(grid[2]))
    widths_ax3 = diff(get_extended_ticks(grid[3]))
    Δw_max_factor = T(1e-2)
    #= 
        Δw_max_factor is chosen by trying out different values for it. 
        This value seems to be okay if the grid is not to unevenly spaced. 
        But this can be secured via the `max_ratio`- and the `the max_tick_distance`-keywords.
        The critical parameter is `csgtol` inside the `in`-method, which is currently defined through 
        the widths of the voxel of the grid point next to the calculated intersection point.
        It would be probably better to turn this into an `NTuple{3,T}` and pass the widths of the voxel instead
        of one single value (`csgtol`). However, inside the `in` methods this would have to be transformed 
        due to rotations and maybe the shape of the primitive. This might not be straight forward.
        But this can be investigated and certainly improved in the future.
    =#
    for i2 in eachindex(ticks[2])
        for i1 in eachindex(ticks[1])
            l = ConstructiveSolidGeometry.Line(CartesianPoint(CylindricalPoint{T}(ticks[1][i1], ticks[2][i2], zero(T))), eZ)
            pts = ConstructiveSolidGeometry.intersection(face, l)
            for pt in pts
                i3 = searchsortednearest(ticks[3], pt[3])
                csgtol = abs(ticks[3][i3] - pt[3])
                Δw_max = Δw_max_factor * max(
                    widths_ax1[i1], widths_ax1[i1+1],
                    # widths_ax2[i2], widths_ax2[i2+1],
                )
                if csgtol < Δw_max
                    csgtol = Δw_max
                end
                if in(pt, geometry, csgtol = csgtol)
                    point_types[i1, i2, i3] = zero(PointType)
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
    #         i2 = searchsortednearest(ticks[2], pt[2])
    #         if in(pt, geometry, csgtol = abs(ticks[2][i2] - pt[2]))
    #             point_types[i1, i2, i3] = zero(PointType)
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
                Δw_max = Δw_max_factor * max(
                    widths_ax3[i3], widths_ax3[i3+1],
                    # widths_ax2[i2], widths_ax2[i2+1],
                )
                if csgtol < Δw_max
                    csgtol = Δw_max
                end
                if abs(pt_cyl[2] - ticks[2][i2]) < 0.1 && in(pt_car, geometry, csgtol = csgtol)
                    point_types[i1, i2, i3] = zero(PointType)
                    potential[i1, i2, i3] = pot_value
                end
            end
        end
    end 
    nothing
end

function mark_bulk_bits!(point_types::Array{PointType, 3})
    i1max, i2max, i3max = size(point_types)
    for i in findall(point_types .& pn_junction_bit .> 0)
        i1, i2, i3 = i[1], i[2], i[3]
        point_types[i1,i2,i3] += bulk_bit * (all([point_types[j1,j2,j3] & pn_junction_bit > 0 
            for j1 in max(i1-1,1):min(i1+1,i1max), j2 in max(i2-1,1):min(i2+1,i2max), j3 in max(i3-1,1):min(i3+1,i3max)]))
    end
    point_types
end