@recipe function f( pss::PotentialSimulationSetup{T, 3, Cylindrical};
                    r = missing,
                    φ = missing,
                    z = missing,
                    n_points_in_φ = 36 ) where {T}
    grid::Grid{T, 3, Cylindrical} = pss.grid
    layout --> (2, 2)

    cross_section::Symbol, idx::Int = if ismissing(φ) && ismissing(r) && ismissing(z)
        :φ, 1
    elseif !ismissing(φ) && ismissing(r) && ismissing(z)
        φ_rad::T = T(deg2rad(φ))
        while !(grid.φ.interval.left <= φ_rad <= grid.φ.interval.right)
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
        grid.φ[idx]
    elseif cross_section == :r
        grid.r[idx]
    elseif cross_section == :z
        grid.z[idx]
    end

    if cross_section == :φ
        φ --> value
    elseif cross_section == :z
        z --> value
    elseif cross_section == :r
        r --> value
    end

    @series begin
        subplot := 1
        PointTypes(pss, n_points_in_φ=n_points_in_φ)
    end
    @series begin
        subplot := 2
        EffectiveChargeDensity(pss)
    end
    @series begin
        subplot := 3
        DielectricDistribution(pss)
    end
    @series begin
        subplot := 4
        ElectricPotential(pss, n_points_in_φ=n_points_in_φ)
    end
end


@recipe function f( pss::PotentialSimulationSetup{T, 3, Cartesian};
                    x = missing,
                    y = missing,
                    z = missing ) where {T}
    grid::Grid{T, 3, Cartesian} = pss.grid
    layout --> (2, 2)

    size --> (1000, 1000)

    cross_section::Symbol, idx::Int = if ismissing(x) && ismissing(y) && ismissing(z)
        :x, searchsortednearest(grid[:x], T(0))
    elseif !ismissing(x) && ismissing(y) && ismissing(z)
        :x, searchsortednearest(grid[:x], T(x))
    elseif ismissing(x) && !ismissing(y) && ismissing(z)
        :y, searchsortednearest(grid[:y], T(y))
    elseif ismissing(x) && ismissing(y) && !ismissing(z)
        :z, searchsortednearest(grid.z, T(z))
    else
        error(ArgumentError, ": Only one of the keywords `x, y, z` is allowed.")
    end

    if cross_section == :x
        x --> grid[:x][idx]
    elseif cross_section == :y
        y --> grid[:y][idx]
    elseif cross_section == :z
        z --> grid.z[idx]
    end


    @series begin
        subplot := 1
        PointTypes(pss)
    end
    @series begin
        subplot := 2
        EffectiveChargeDensity(pss)
    end
    @series begin
        subplot := 3
        DielectricDistribution(pss)
    end
    @series begin
        subplot := 4
        ElectricPotential(pss)
    end
end
