function _get_potential_plot_information(p::ElectricPotential)::Tuple{Symbol, Tuple{Quantity, Quantity}, String, Unitful.Units}
    return (:viridis, extrema(p.data).*internal_voltage_unit, "Electric Potential", internal_voltage_unit)
end

function _get_potential_plot_information(::WeightingPotential{T})::Tuple{Symbol, Tuple{T,T}, String, Unitful.Units} where {T}
    return (:viridis, (zero(T), one(T)), "Weighting Potential", Unitful.NoUnits)
end

function _get_potential_plot_information(ecd::EffectiveChargeDensity)::Tuple{Symbol, Tuple{Quantity, Quantity}, String, Unitful.Units}
    return (:inferno, extrema(ecd.data).*u"V/m^2", "Effective Charge Density", u"V/m^2")
end

function _get_potential_plot_information(::PointTypes{T})::Tuple{Symbol, Tuple{T,T}, String, Unitful.Units} where {T}
    return (:viridis, (zero(0), 7*one(T)), "Point Type Map", Unitful.NoUnits)
end

function _get_potential_plot_information(::ImpurityScale{T})::Tuple{Symbol, Tuple{T,T}, String, Unitful.Units} where {T}
    return (:viridis, (zero(T), one(T)), "Impurity Scale", Unitful.NoUnits)
end

function get_crosssection_idx_and_value(grid::Grid, ::Any, ::Any, ::Any)
    error(ArgumentError, ": The cross section your are trying to plot is not given in the right format.\nPlease pass a number with matching or no units.")
end


## Potentials with cylindrical grids

function get_crosssection_idx_and_value(grid::CylindricalGrid{T}, r::Union{Real, LengthQuantity, Missing}, φ::Union{Real, AngleQuantity, Missing}, z::Union{Real, LengthQuantity, Missing})::Tuple{Symbol,Int,T,Unitful.Units} where {T <: SSDFloat}

    cross_section::Symbol, idx::Int, units::Unitful.Units = if ismissing(φ) && ismissing(r) && ismissing(z)
        return get_crosssection_idx_and_value(grid, r, T(0.0), z)
    elseif !ismissing(φ) && ismissing(r) && ismissing(z)
        φ_rad::T = T(to_internal_units(φ isa Real ? deg2rad(φ) : φ))
        while !(grid.φ.interval.left <= φ_rad <= grid.φ.interval.right) && grid.φ.interval.right != grid.φ.interval.left
            if φ_rad > grid.φ.interval.right
                φ_rad -= width(grid.φ.interval)
            elseif φ_rad < grid.φ.interval.left
                φ_rad += width(grid.φ.interval)
            end
        end
        :φ, searchsortednearest(grid.φ, φ_rad), φ isa Real ? u"°" : unit(φ)
    elseif ismissing(φ) && !ismissing(r) && ismissing(z)
        :r, searchsortednearest(grid.r, T(to_internal_units(r))), r isa Real ? internal_length_unit : unit(r)
    elseif ismissing(φ) && ismissing(r) && !ismissing(z)
        :z, searchsortednearest(grid.z, T(to_internal_units(z))), z isa Real ? internal_length_unit : unit(z)
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

    cross_section, idx, value, units
end

@recipe function f(sp::ScalarPotential{T,3,Cylindrical}; r = missing, φ = missing, z = missing, contours_equal_potential = false, full_det = false) where {T <: SSDFloat}

    gradient::Symbol, clims::Tuple{RealQuantity, RealQuantity}, name::String, punit::Unitful.Units = _get_potential_plot_information(sp)

    if !(sp.grid[2][end] - sp.grid[2][1] ≈ 2π) sp = get_2π_potential(sp, n_points_in_φ = 72) end

    grid::CylindricalGrid{T} = sp.grid
    cross_section::Symbol, idx::Int, value::T, units::Unitful.Units = get_crosssection_idx_and_value(grid, r, φ, z)
    
    title --> name * " @ $(cross_section) = $(round(units, Float64(uconvert(units,value*(cross_section == :φ ? u"°" : internal_length_unit))), sigdigits=3))"
    colorbar_title --> name
    seriescolor --> gradient
    clims --> clims
    unitformat --> :slash

    sp, cross_section, idx, value, punit, contours_equal_potential, full_det
end

@recipe function f(sp::ScalarPotential{T,3,Cylindrical}, cross_section::Symbol, idx::Int, value::T, punit::Unitful.Units, contours_equal_potential::Bool, full_det::Bool) where {T <: SSDFloat}

    grid::CylindricalGrid{T} = sp.grid
    data = sp.data
    if eltype(data) == PointType
        data = data .& 0x07
    end
    @series begin
        seriestype := :heatmap
        foreground_color_border --> nothing
        tick_direction --> :out
        unitformat --> :slash
        if cross_section == :φ
            aspect_ratio --> 1
            xguide --> "r"
            xunit --> internal_length_unit
            xlims --> extrema(grid.r).*internal_length_unit
            yguide --> "z"
            yunit --> internal_length_unit
            ylims --> extrema(grid.z).*internal_length_unit
            gr_ext::Array{T,1} = midpoints(get_extended_ticks(grid.r))
            gz_ext::Array{T,1} = midpoints(get_extended_ticks(grid.z))
            if full_det
                cross_section_dummy, idx_mirror, value_dummy = get_crosssection_idx_and_value(grid, missing, value+180, missing)
                extended_data =  cat(data[end:-1:2, idx_mirror, :]', data[:, idx, :]', dims = 2)
                xlims := (-1*grid.r[end],grid.r[end]).*internal_length_unit
                vcat(-1 .* grid.r[end:-1:2], grid.r)*internal_length_unit, grid.z*internal_length_unit, extended_data*punit
             else
                midpoints(gr_ext)*internal_length_unit, midpoints(gz_ext)*internal_length_unit, data[:,idx,:]'*punit
            end
        elseif cross_section == :r
            xguide --> "φ"
            xunit --> internal_angle_unit
            yguide --> "z"
            yunit --> internal_length_unit
            ylims --> extrema(grid.z)*internal_length_unit
            grid.φ*internal_angle_unit, grid.z*internal_length_unit, data[idx,:,:]'*punit
        elseif cross_section == :z
            projection --> :polar
            xguide --> ""
            yguide --> ""
            grid.φ*internal_angle_unit, grid.r*internal_length_unit, data[:,:,idx].*punit
        end
    end

    if contours_equal_potential && cross_section == :φ
        @series begin
            seriescolor := :thermal
            seriestype := :contours
            unitformat --> :slash
            aspect_ratio --> 1
            xguide := "r"
            yguide := "z"
            xlims --> extrema(grid.r).*internal_length_unit
            ylims --> extrema(grid.z).*internal_length_unit
            if full_det
                cross_section_dummy, idx_mirror, value_dummy = get_crosssection_idx_and_value(grid, missing, value+180, missing)
                extended_data =  cat(data[end:-1:2, idx_mirror, :]', data[:, idx, :]', dims = 2)
                xlims := (-last(grid.r),last(grid.r)).*internal_length_unit
                vcat(-1 .* grid.r[end:-1:2], grid.r)*internal_length_unit, grid.z*internal_length_unit, extended_data.*punit
             else
                # midpoints(gr_ext), midpoints(gz_ext), data[:,idx,:]'*unit
                grid.r*internal_length_unit, grid.z*internal_length_unit, data[:,idx,:]'.*punit
            end
        end
    end
end



@recipe function f(ϵ::DielectricDistribution{T,3,Cylindrical}; r = missing, φ = missing, z = missing) where {T <: SSDFloat}

    grid::CylindricalGrid{T} = ϵ.grid
    cross_section::Symbol, idx::Int, value::T, units::Unitful.Units = get_crosssection_idx_and_value(grid, r, φ, z)

    seriestype --> :heatmap
    seriescolor --> :inferno
    foreground_color_border --> nothing
    tick_direction --> :out
    title --> "Dielectric Distribution @ $(cross_section) = $(round(units, Float64(uconvert(units,value*(cross_section == :φ ? u"°" : internal_length_unit))), sigdigits=3))"
    colorbar_title --> "Dielectric Distribution"
    unitformat --> :slash

    gr_ext::Array{T,1} = midpoints(get_extended_ticks(grid.r))
    gφ_ext::Array{T,1} = midpoints(get_extended_ticks(grid.φ))
    gz_ext::Array{T,1} = midpoints(get_extended_ticks(grid.z))

    @series begin
        if cross_section == :φ
            aspect_ratio --> 1
            xguide --> "r"
            xunit --> internal_length_unit
            yguide --> "z"
            yunit --> internal_length_unit
            xlims --> (grid.r[2],grid.r[end-1]).*internal_length_unit
            ylims --> (grid.z[2],grid.z[end-1]).*internal_length_unit
            gr_ext*internal_length_unit, gz_ext*internal_length_unit, ϵ.data[:,idx,:]'
        elseif cross_section == :r
            xguide --> "φ"
            xunit --> internal_angle_unit
            yguide --> "z"
            yunit --> internal_length_unit
            ylims --> (grid.z[2],grid.z[end-1]).*internal_length_unit
            gφ_ext*internal_angle_unit, gz_ext*internal_length_unit, ϵ.data[idx,:,:]'
        elseif cross_section == :z
            projection --> :polar
            gφ_ext*internal_angle_unit, gr_ext[2:end-1]*internal_length_unit, ϵ.data[2:end-1,:,idx]
        end
    end
end



## Potentials with Cartesian grids

function get_crosssection_idx_and_value(grid::CartesianGrid3D{T}, x::Union{Real, LengthQuantity, Missing}, y::Union{Real, LengthQuantity, Missing}, z::Union{Real, LengthQuantity, Missing})::Tuple{Symbol,Int,T,Unitful.Units} where {T <: SSDFloat}

    cross_section::Symbol, idx::Int, units::Unitful.Units = if ismissing(x) && ismissing(y) && ismissing(z)
        return get_crosssection_idx_and_value(grid, T(0.0), y, z)
    elseif !ismissing(x) && ismissing(y) && ismissing(z)
        :x, searchsortednearest(grid.x, T(to_internal_units(x))), x isa Real ? internal_length_unit : unit(x)
    elseif ismissing(x) && !ismissing(y) && ismissing(z)
        :y, searchsortednearest(grid.y, T(to_internal_units(y))), y isa Real ? internal_length_unit : unit(y)
    elseif ismissing(x) && ismissing(y) && !ismissing(z)
        :z, searchsortednearest(grid.z, T(to_internal_units(z))), z isa Real ? internal_length_unit : unit(z)
    else
        error(ArgumentError, ": Only one of the keywords `x, y, z` is allowed.")
    end
    value::T = grid[cross_section][idx]
    cross_section, idx, value, units
end

@recipe function f(sp::ScalarPotential{T,3,Cartesian}; x = missing, y = missing, z = missing, contours_equal_potential = false) where {T <: SSDFloat}
    
    gradient::Symbol, clims::Tuple{RealQuantity, RealQuantity}, name::String, punit::Unitful.Units = _get_potential_plot_information(sp)

    grid::CartesianGrid3D{T} = sp.grid
    cross_section::Symbol, idx::Int, value::T, units::Unitful.Units = get_crosssection_idx_and_value(grid, x, y, z)
    
    title --> name * " @ $(cross_section) = $(round(units, Float64(uconvert(units,value*internal_length_unit)), sigdigits=3))"
    colorbar_title --> name
    seriescolor --> gradient
    clims --> clims

    sp, cross_section, idx, value, punit, contours_equal_potential
end

@recipe function f(sp::ScalarPotential{T,3,Cartesian}, cross_section::Symbol, idx::Int, value::T, punit::Unitful.Units, contours_equal_potential::Bool) where {T <: SSDFloat}

    grid::CartesianGrid3D{T} = sp.grid
    data = sp.data
    if eltype(data) == PointType
        data = data .& 0x07
    end
    @series begin
        seriestype := :heatmap
        foreground_color_border --> nothing
        tick_direction --> :out
        unitformat --> :slash
        xunit --> internal_length_unit
        yunit --> internal_length_unit
        if cross_section == :x
            aspect_ratio --> 1
            xguide --> "y"
            yguide --> "z"
            xlims --> extrema(grid.y).*internal_length_unit
            ylims --> extrema(grid.z).*internal_length_unit
            gy_ext = midpoints(get_extended_ticks(grid.y))
            gz_ext = midpoints(get_extended_ticks(grid.z))
            midpoints(gy_ext)*internal_length_unit, midpoints(gz_ext)*internal_length_unit, data[idx,:,:]'*punit
        elseif cross_section == :y
            aspect_ratio --> 1
            xguide --> "x"
            yguide --> "z"
            xlims --> extrema(grid.x).*internal_length_unit
            ylims --> extrema(grid.z).*internal_length_unit
            gx_ext = midpoints(get_extended_ticks(grid.x))
            gz_ext = midpoints(get_extended_ticks(grid.z))
            midpoints(gx_ext)*internal_length_unit, midpoints(gz_ext)*internal_length_unit, data[:,idx,:]'*punit
        elseif cross_section == :z
            aspect_ratio --> 1
            xguide --> "x"
            yguide --> "y"
            xlims --> extrema(grid.x).*internal_length_unit
            ylims --> extrema(grid.y).*internal_length_unit
            gx_ext = midpoints(get_extended_ticks(grid.x))
            gy_ext = midpoints(get_extended_ticks(grid.y))
            midpoints(gx_ext)*internal_length_unit, midpoints(gy_ext)*internal_length_unit, data[:,:,idx]'*punit
        end
    end

    if contours_equal_potential
        @series begin
            seriescolor := :thermal
            seriestype := :contours
            aspect_ratio --> 1
            unitformat --> :slash
            if cross_section == :x
                gy_ext = midpoints(get_extended_ticks(grid.y))
                gz_ext = midpoints(get_extended_ticks(grid.z))
                midpoints(gy_ext)*internal_length_unit, midpoints(gz_ext)*internal_length_unit, data[idx,:,:]'*punit
            elseif cross_section == :y
                gx_ext = midpoints(get_extended_ticks(grid.x))
                gz_ext = midpoints(get_extended_ticks(grid.z))
                midpoints(gx_ext)*internal_length_unit, midpoints(gz_ext)*internal_length_unit, data[:,idx,:]'*punit
            elseif cross_section == :z
                gx_ext = midpoints(get_extended_ticks(grid.x))
                gy_ext = midpoints(get_extended_ticks(grid.y))
                midpoints(gx_ext)*internal_length_unit, midpoints(gy_ext)*internal_length_unit, data[:,:,idx]'*punit
            end
        end
    end
end


@recipe function f(ϵ::DielectricDistribution{T,3,Cartesian}; x = missing, y = missing, z = missing) where {T <: SSDFloat}

    grid::CartesianGrid3D{T} = ϵ.grid
    cross_section::Symbol, idx::Int, value::T, units::Unitful.Units = get_crosssection_idx_and_value(grid, x, y, z)
    seriestype --> :heatmap
    seriescolor --> :inferno
    foreground_color_border --> nothing
    tick_direction --> :out
    title --> "Dielectric Distribution @ $(cross_section) = $(round(units, Float64(uconvert(units,value*internal_length_unit)), sigdigits=3))"
    colorbar_title --> "Dielectric Distribution"

    gx_ext::Array{T,1} = midpoints(get_extended_ticks(grid.x))
    gy_ext::Array{T,1} = midpoints(get_extended_ticks(grid.y))
    gz_ext::Array{T,1} = midpoints(get_extended_ticks(grid.z))

    @series begin
        unitformat --> :slash
        xunit --> internal_length_unit
        yunit --> internal_length_unit
        if cross_section == :x
            aspect_ratio --> 1
            xguide --> "y"
            yguide --> "z"
            xlims --> extrema(grid.y).*internal_length_unit
            ylims --> extrema(grid.z).*internal_length_unit
            gy_ext*internal_length_unit, gz_ext*internal_length_unit, ϵ.data[idx,:,:]'
        elseif cross_section == :y
            aspect_ratio --> 1
            xguide --> "x"
            yguide --> "z"
            xlims --> extrema(grid.x).*internal_length_unit
            ylims --> extrema(grid.z).*internal_length_unit
            gx_ext*internal_length_unit, gz_ext*internal_length_unit, ϵ.data[:,idx,:]'
        elseif cross_section == :z
            aspect_ratio --> 1
            xguide --> "x"
            yguide --> "y"
            xlims --> extrema(grid.x).*internal_length_unit
            ylims --> extrema(grid.y).*internal_length_unit
            gx_ext*internal_length_unit, gy_ext*internal_length_unit, ϵ.data[:,:,idx]'
        end
    end
end
