function get_crosssection_idx_and_value(grid::CylindricalGrid{T}, r, φ, z)::Tuple{Symbol,Int,T} where {T <: SSDFloat}

    cross_section::Symbol, idx::Int = if ismissing(φ) && ismissing(r) && ismissing(z)
        return get_crosssection_idx_and_value(grid, r, T(0.0), z)
    elseif !ismissing(φ) && ismissing(r) && ismissing(z)
        φ_rad::T = T(deg2rad(φ))
        while !(grid.φ.interval.left <= φ_rad <= grid.φ.interval.right) && grid.φ.interval.right != grid.φ.interval.left
            if φ_rad > grid.φ.interval.right
                φ_rad -= width(grid.φ.interval)
            elseif φ_rad < grid.φ.interval.left
                φ_rad += width(grid.φ.interval)
            end
        end
        :φ, searchsortednearest(grid.φ, φ_rad)
    elseif ismissing(φ) && !ismissing(r) && ismissing(z)
        :r, searchsortednearest(grid.r, T(r))
    elseif ismissing(φ) && ismissing(r) && !ismissing(z)
        :z, searchsortednearest(grid.z, T(z))
    else
        error(ArgumentError, ": Only one of the keywords `r, φ, z` is allowed.")
    end

    value::T = if cross_section == :φ
        rad2deg(grid.φ[idx])
    elseif cross_section == :r
        grid.r[idx]
    elseif cross_section == :z
        grid.z[idx]
    end

    cross_section, idx, value
end



@recipe function f(epot::ElectricPotential{T,3,Cylindrical}; r = missing, φ = missing, z = missing, contours_equal_potential = false, full_det = false) where {T <: SSDFloat}

    if !(epot.grid[2][end] - epot.grid[2][1] ≈ 2π) epot = get_2π_potential(epot, n_points_in_φ = 72) end

    grid::CylindricalGrid{T} = epot.grid
    cross_section::Symbol, idx::Int, value::T = get_crosssection_idx_and_value(grid, r, φ, z)

    seriescolor --> :viridis
    title --> "Electric Potential @ $(cross_section) = $(round(value,sigdigits=2))"*(cross_section == :φ ? "°" : "m")
    colorbar_title --> "ElectricPotential"

    epot, cross_section, idx, value, internal_voltage_unit, contours_equal_potential, full_det
end


@recipe function f(wpot::WeightingPotential{T,3,Cylindrical}; r = missing, φ = missing, z = missing, contours_equal_potential = false, full_det = false) where {T <: SSDFloat}

    if !(wpot.grid[2][end] - wpot.grid[2][1] ≈ 2π)
        wpot = get_2π_potential(wpot, n_points_in_φ = 72)
    end

    grid::CylindricalGrid{T} = wpot.grid
    cross_section::Symbol, idx::Int, value::T = get_crosssection_idx_and_value(grid, r, φ, z)

    seriescolor --> :viridis
    clims --> (0,1)
    title --> "Weighting Potential @ $(cross_section) = $(round(value,sigdigits=2))"*(cross_section == :φ ? "°" : "m")

    wpot, cross_section, idx, value, Unitful.NoUnits, contours_equal_potential, full_det
end


@recipe function f(ρ::EffectiveChargeDensity{T,3,Cylindrical}; r = missing, φ = missing, z = missing, full_det = false) where {T <: SSDFloat}

    if !(ρ.grid[2][end] - ρ.grid[2][1] ≈ 2π)
        ρ = get_2π_potential(ρ, n_points_in_φ = 72)
    end

    grid::CylindricalGrid{T} = ρ.grid
    cross_section::Symbol, idx::Int, value::T = get_crosssection_idx_and_value(grid, r, φ, z)

    seriescolor --> :inferno
    title --> "Effective Charge Density @ $(cross_section) = $(round(value,sigdigits=2))"*(cross_section == :φ ? "°" : "m")
    ρ, cross_section, idx, value, u"V*m", false, full_det
end


@recipe function f(point_types::PointTypes{T,3,Cylindrical}; r = missing, φ = missing, z = missing, full_det = false) where {T <: SSDFloat}

    if !(point_types.grid[2][end] - point_types.grid[2][1] ≈ 2π)
        point_types = get_2π_potential(point_types, n_points_in_φ = 72)
    end

    grid::CylindricalGrid{T} = point_types.grid
    cross_section::Symbol, idx::Int, value::T = get_crosssection_idx_and_value(grid, r, φ, z)

    seriescolor --> :viridis
    clims --> (0,7)
    title --> "Point Type Map @ $(cross_section) = $(round(value,sigdigits=2))"*(cross_section == :φ ? "°" : "m")

    point_types, cross_section, idx, value, Unitful.NoUnits, false, full_det
end

@recipe function f(sp::ScalarPotential{T,3,Cylindrical}, cross_section::Symbol, idx::Int, value::T, unit::Unitful.Units, contours_equal_potential::Bool = false, full_det::Bool = false) where {T <: SSDFloat}
    grid::CylindricalGrid{T} = sp.grid
    @series begin
        seriestype := :heatmap
        foreground_color_border --> nothing
        tick_direction --> :out
        if cross_section == :φ
            aspect_ratio --> 1
            xguide --> "r"
            yguide --> "z"
            unitformat --> :slash
            xlims --> (grid.r[1],grid.r[end])
            ylims --> (grid.z[1],grid.z[end])
            gr_ext::Array{T,1} = midpoints(get_extended_ticks(grid.r))
            gz_ext::Array{T,1} = midpoints(get_extended_ticks(grid.z))
            if full_det
                cross_section_dummy, idx_mirror, value_dummy = get_crosssection_idx_and_value(grid, missing, value+180, missing)
                extended_data =  cat(sp.data[end:-1:2, idx_mirror, :]', sp.data[:, idx, :]', dims = 2)
                xlims := (-1*grid.r[end],grid.r[end])
                vcat(-1 .* grid.r[end:-1:2], grid.r)*internal_length_unit, grid.z*internal_length_unit, extended_data*unit
             else
                midpoints(gr_ext)*internal_length_unit, midpoints(gz_ext)*internal_length_unit, sp.data[:,idx,:]'*unit
            end
        elseif cross_section == :r
            xguide --> "φ"
            yguide --> "z"
            ylims --> (grid.z[1],grid.z[end])
            grid.φ*internal_angle_unit, grid.z*internal_length_unit, sp.data[idx,:,:]'*unit
        elseif cross_section == :z
            projection --> :polar
            xguide --> ""
            yguide --> ""
            grid.φ*internal_angle_unit, grid.r*internal_length_unit, sp.data[:,:,idx]*unit
        end
    end

    if contours_equal_potential && cross_section == :φ
        @series begin
            seriescolor := :thermal
            seriestype := :contours
            #if cross_section == :φ
                aspect_ratio --> 1
                xguide := "r"
                yguide := "z"
                unitformat --> :slash
                xlims --> (grid.r[1],grid.r[end])
                ylims --> (grid.z[1],grid.z[end])
                if full_det
                    cross_section_dummy, idx_mirror, value_dummy = get_crosssection_idx_and_value(grid, missing, value+180, missing)
                    extended_data =  cat(sp.data[end:-1:2, idx_mirror, :]', sp.data[:, idx, :]', dims = 2)
                    xlims := (-1*grid.r[end],grid.r[end])
                    vcat(-1 .* grid.r[end:-1:2], grid.r)*internal_length_unit, grid.z*internal_length_unit, extended_data*unit
                 else
                    # midpoints(gr_ext), midpoints(gz_ext), sp.data[:,idx,:]'*unit
                    grid.r*internal_length_unit, grid.z*internal_length_unit, sp.data[:,idx,:]'*unit
                end
                #=
            elseif cross_section == :r
                xguide --> "φ / °"
                yguide --> "z / m"
                ylims --> (grid.z[1],grid.z[end])
                grid.φ, grid.z, sp.data[idx,:,:]'
            elseif cross_section == :z
                projection --> :polar
                grid.φ, grid.r, sp.data[:,:,idx]
            end
            =#
        end
    end
end



@recipe function f(ϵ::DielectricDistribution{T,3,Cylindrical}; φ = 0) where {T <: SSDFloat}

    grid::CylindricalGrid{T} = ϵ.grid
    cross_section::Symbol, idx::Int, value::T = get_crosssection_idx_and_value(grid, missing, φ, missing)

    seriestype --> :heatmap
    seriescolor --> :inferno
    foreground_color_border --> nothing
    tick_direction --> :out
    title --> "Dielectric Distribution @ $(cross_section) = $(round(value,sigdigits=2))"*(cross_section == :φ ? "°" : "m")

    gr_ext::Array{T,1} = midpoints(get_extended_ticks(grid.r))
    gφ_ext::Array{T,1} = midpoints(get_extended_ticks(grid.φ))
    gz_ext::Array{T,1} = midpoints(get_extended_ticks(grid.z))

    @series begin
        if cross_section == :φ
            aspect_ratio --> 1
            xguide --> "r"
            yguide --> "z"
            unitformat --> :slash
            xlims --> (grid.r[1],grid.r[end])
            ylims --> (grid.z[1],grid.z[end])
            gr_ext*internal_length_unit, gz_ext*internal_length_unit, ϵ.data[:,idx,:]'
        elseif cross_section == :r
            xguide --> "φ"
            yguide --> "z"
            ylims --> (grid.z[1],grid.z[end])
            rad2deg_backend.(gφ_ext)*internal_angle_unit, gz_ext*internal_length_unit, ϵ.data[idx,:,:]'
        elseif cross_section == :z
            projection --> :polar
            rad2deg_backend.(gφ_ext)*internal_angle_unit, gr_ext*internal_length_unit, ϵ.data[:,:,idx]
        end
    end
end



function get_crosssection_idx_and_value(grid::CartesianGrid3D{T}, x, y, z)::Tuple{Symbol,Int,T} where {T <: SSDFloat}

    cross_section::Symbol, idx::Int = if ismissing(x) && ismissing(y) && ismissing(z)
        return get_crosssection_idx_and_value(grid, T(0.0), y, z)
    elseif !ismissing(x) && ismissing(y) && ismissing(z)
        :x, searchsortednearest(grid.x, T(x))
    elseif ismissing(x) && !ismissing(y) && ismissing(z)
        :y, searchsortednearest(grid.y, T(y))
    elseif ismissing(x) && ismissing(y) && !ismissing(z)
        :z, searchsortednearest(grid.z, T(z))
    else
        error(ArgumentError, ": Only one of the keywords `x, y, z` is allowed.")
    end
    value::T = grid[cross_section][idx]
    cross_section, idx, value
end

@recipe function f(epot::ElectricPotential{T,3,Cartesian}; x = missing, y = missing, z = missing, contours_equal_potential = false) where {T <: SSDFloat}
    grid::CartesianGrid3D{T} = epot.grid
    cross_section::Symbol, idx::Int, value::T = get_crosssection_idx_and_value(grid, x, y, z)

    seriescolor --> :viridis
    title --> "Electric Potential @ $(cross_section) = $(round(value,sigdigits=2))m"

    epot, cross_section, idx, value, internal_voltage_unit, contours_equal_potential
end

@recipe function f(wpot::WeightingPotential{T,3,Cartesian}; x = missing, y = missing, z = missing, contours_equal_potential = false) where {T <: SSDFloat}

    grid::CartesianGrid3D{T} = wpot.grid
    cross_section::Symbol, idx::Int, value::T = get_crosssection_idx_and_value(grid, x, y, z)

    seriescolor --> :viridis
    clims --> (0,1)
    title --> "Weighting Potential @ $(cross_section) = $(round(value,sigdigits=2))m"

    wpot, cross_section, idx, value, Unitful.NoUnits, contours_equal_potential
end


@recipe function f(ρ::EffectiveChargeDensity{T,3,Cartesian}; x = missing, y = missing, z = missing) where {T <: SSDFloat}
    grid::CartesianGrid3D{T} = ρ.grid
    cross_section::Symbol, idx::Int, value::T = get_crosssection_idx_and_value(grid, x, y, z)

    seriescolor --> :inferno
    title --> "Effective Charge Density @ $(cross_section) = $(round(value,sigdigits=2))m"

    ρ, cross_section, idx, value, u"V*m"
end


@recipe function f(point_types::PointTypes{T,3,Cartesian}; x = missing, y = missing, z = missing) where {T <: SSDFloat}

    grid::CartesianGrid3D{T} = point_types.grid
    cross_section::Symbol, idx::Int, value::T = get_crosssection_idx_and_value(grid, x, y, z)

    seriescolor --> :viridis
    clims --> (0,7)
    title --> "Point Type Map @ $(cross_section) = $(round(value,sigdigits=2))m"

    point_types, cross_section, idx, value, Unitful.NoUnits
end


@recipe function f(sp::ScalarPotential{T,3,Cartesian}, cross_section::Symbol, idx::Int, value::T, unit::Unitful.Units, contours_equal_potential::Bool = false) where {T <: SSDFloat}
    grid::CartesianGrid3D{T} = sp.grid
    @series begin
        seriestype := :heatmap
        foreground_color_border --> nothing
        tick_direction --> :out
        unitformat --> :slash
        if cross_section == :x
            aspect_ratio --> 1
            xguide --> "y"
            yguide --> "z"
            xlims --> (grid.y[1],grid.y[end])
            ylims --> (grid.z[1],grid.z[end])
            gy_ext = midpoints(get_extended_ticks(grid.y))
            gz_ext = midpoints(get_extended_ticks(grid.z))
            midpoints(gy_ext)*internal_length_unit, midpoints(gz_ext)*internal_length_unit, sp.data[idx,:,:]'*unit
        elseif cross_section == :y
            aspect_ratio --> 1
            xguide --> "x"
            yguide --> "z"
            xlims --> (grid.x[1],grid.x[end])
            ylims --> (grid.z[1],grid.z[end])
            gx_ext = midpoints(get_extended_ticks(grid.x))
            gz_ext = midpoints(get_extended_ticks(grid.z))
            midpoints(gx_ext)*internal_length_unit, midpoints(gz_ext)*internal_length_unit, sp.data[:,idx,:]'*unit
        elseif cross_section == :z
            aspect_ratio --> 1
            xguide --> "x"
            yguide --> "y"
            xlims --> (grid.x[1],grid.x[end])
            ylims --> (grid.y[1],grid.y[end])
            gx_ext = midpoints(get_extended_ticks(grid.x))
            gy_ext = midpoints(get_extended_ticks(grid.y))
            midpoints(gx_ext)*internal_length_unit, midpoints(gy_ext)*internal_length_unit, sp.data[:,:,idx]'*unit
        end
    end

    if contours_equal_potential
        @series begin
            seriescolor := :thermal
            seriestype := :contours
            aspect_ratio --> 1
            unitformat --> :slash
            if cross_section == :x
                xguide --> "y"
                yguide --> "z"
                xlims --> (grid.y[1],grid.y[end])
                ylims --> (grid.z[1],grid.z[end])
                gy_ext = midpoints(get_extended_ticks(grid.y))
                gz_ext = midpoints(get_extended_ticks(grid.z))
                midpoints(gy_ext)*internal_length_unit, midpoints(gz_ext)*internal_length_unit, sp.data[idx,:,:]'*unit
            elseif cross_section == :y
                xguide --> "x"
                yguide --> "z"
                xlims --> (grid.x[1],grid.x[end])
                ylims --> (grid.z[1],grid.z[end])
                gx_ext = midpoints(get_extended_ticks(grid.x))
                gz_ext = midpoints(get_extended_ticks(grid.z))
                midpoints(gx_ext)*internal_length_unit, midpoints(gz_ext)*internal_length_unit, sp.data[:,idx,:]'*unit
            elseif cross_section == :z
                xguide --> "x"
                yguide --> "y"
                xlims --> (grid.x[1],grid.x[end])
                ylims --> (grid.y[1],grid.y[end])
                gx_ext = midpoints(get_extended_ticks(grid.x))
                gy_ext = midpoints(get_extended_ticks(grid.y))
                midpoints(gx_ext)*internal_length_unit, midpoints(gy_ext)*internal_length_unit, sp.data[:,:,idx]'*unit
            end
        end
    end
end


@recipe function f(ϵ::DielectricDistribution{T,3,Cartesian}; x = missing, y = missing, z = missing) where {T <: SSDFloat}

    grid::CartesianGrid3D{T} = ϵ.grid
    cross_section::Symbol, idx::Int, value::T = get_crosssection_idx_and_value(grid, x, y, z)

    seriestype --> :heatmap
    seriescolor --> :inferno
    foreground_color_border --> nothing
    tick_direction --> :out
    title --> "Dielectric Distribution @ $(cross_section) = $(round(value,sigdigits=2))m"

    gx_ext::Array{T,1} = midpoints(get_extended_ticks(grid.x))
    gy_ext::Array{T,1} = midpoints(get_extended_ticks(grid.y))
    gz_ext::Array{T,1} = midpoints(get_extended_ticks(grid.z))

    @series begin
        unitformat --> :slash
        if cross_section == :x
            aspect_ratio --> 1
            xguide --> "y"
            yguide --> "z"
            xlims --> (grid.y[1],grid.y[end])
            ylims --> (grid.z[1],grid.z[end])
            gy_ext*internal_length_unit, gz_ext*internal_length_unit, ϵ.data[idx,:,:]'
        elseif cross_section == :y
            aspect_ratio --> 1
            xguide --> "x"
            yguide --> "z"
            xlims --> (grid.x[1],grid.x[end])
            ylims --> (grid.z[1],grid.z[end])
            gx_ext*internal_length_unit, gz_ext*internal_length_unit, ϵ.data[:,idx,:]'
        elseif cross_section == :z
            aspect_ratio --> 1
            xguide --> "x"
            yguide --> "y"
            xlims --> (grid.x[1],grid.x[end])
            ylims --> (grid.y[1],grid.y[end])
            gx_ext*internal_length_unit, gy_ext*internal_length_unit, ϵ.data[:,:,idx]'
        end
    end
end
