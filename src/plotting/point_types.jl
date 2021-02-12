

@recipe function f( pts::PointTypes{T, 3, :cylindrical};
                    r = missing,
                    φ = missing,
                    z = missing ) where {T}
    g::Grid{T, 3, :cylindrical} = pts.grid
   
    seriescolor --> :viridis
    st --> :heatmap
    aspect_ratio --> 1
    
    cross_section::Symbol, idx::Int = if ismissing(φ) && ismissing(r) && ismissing(z)
        :φ, 1
    elseif !ismissing(φ) && ismissing(r) && ismissing(z)
        φ_rad::T = T(deg2rad(φ))
        while !(g.φ.interval.left <= φ_rad <= g.φ.interval.right)
            if φ_rad > g.φ.interval.right
                φ_rad -= g.φ.interval.right - g.φ.interval.left
            elseif φ_rad < g.φ.interval.left
                φ_rad += g.φ.interval.right - g.φ.interval.left
            end
        end
        :φ, searchsortednearest(g.φ, φ_rad)
    elseif ismissing(φ) && !ismissing(r) && ismissing(z)
        :r, searchsortednearest(g.r, T(r))
    elseif ismissing(φ) && ismissing(r) && !ismissing(z)
        :z, searchsortednearest(g.z, T(z))
    else
        error(ArgumentError, ": Only one of the keywords `r, φ, z` is allowed.")
    end
    value::T = if cross_section == :φ
        g.φ[idx]
    elseif cross_section == :r    
        g.r[idx]
    elseif cross_section == :z
        g.z[idx]
    end
    
    @series begin
        clims --> (0, max_pointtype_value)
        if cross_section == :φ
            title --> "Point Type Map @$(cross_section) = $(round(rad2deg(value), sigdigits = 2))"
            xlabel --> "r / m"
            ylabel --> "z / m"
            size --> ( 400, 350 / (g.r[end] - g.r[1]) * (g.z[end] - g.z[1]) )
            g.r, g.z, pts.data[:, idx,:]'
        elseif cross_section == :r
            title --> "Point Type Map @$(cross_section) = $(round(value, sigdigits = 2))"
            g.φ, g.z, pts.data[idx,:,:]'
        elseif cross_section == :z
            title --> "Point Type Map @$(cross_section) = $(round(value, sigdigits = 2))"
            proj --> :polar
            g.φ, g.r, pts.data[:,:,idx]
        end
    end
end



@recipe function f( pts::PointTypes{T, 3, :cartesian};
                    x = missing,
                    y = missing,
                    z = missing ) where {T}
    g::Grid{T, 3, :cartesian} = pts.grid
   
    seriescolor --> :viridis
    st --> :heatmap
    aspect_ratio --> 1
    
    cross_section::Symbol, idx::Int = if ismissing(x) && ismissing(y) && ismissing(z)
        :x, 1
    elseif !ismissing(x) && ismissing(y) && ismissing(z)
        :x, searchsortednearest(g[:x], T(x))
    elseif ismissing(x) && !ismissing(y) && ismissing(z)
        :y, searchsortednearest(g[:y], T(y))
    elseif ismissing(x) && ismissing(y) && !ismissing(z)
        :z, searchsortednearest(g.z, T(z))
    else
        error(ArgumentError, ": Only one of the keywords `r, y, z` is allowed.")
    end
    value::T = if cross_section == :x
        g[:x][idx]
    elseif cross_section == :y
        g[:y][idx]
    elseif cross_section == :z
        g.z[idx]
    end
    
    @series begin
        clims --> (0, max_pointtype_value)
        title --> "Point Type Map @$(cross_section) = $(round(value, sigdigits = 2))"
        if cross_section == :x
            xlabel --> "y / m"
            ylabel --> "z / m"
            g[:y], g.z, pts.data[idx, :, :]'
        elseif cross_section == :y
            xlabel --> "x / m"
            ylabel --> "z / m"
            g[:x], g.z, pts.data[:, idx, :]'
        elseif cross_section == :z
            xlabel --> "x / m"
            ylabel --> "y / m"
            g[:x], g[:y], pts.data[:,:,idx]'
        end
    end
end
