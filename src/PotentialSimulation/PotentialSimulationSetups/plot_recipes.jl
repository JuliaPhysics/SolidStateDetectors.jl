
@recipe function f( pss::PotentialSimulationSetup{T, 3, :cylindrical};
                    r = missing,
                    φ = 0,
                    z = missing,
                    n_points_in_φ = 36 ) where {T}
    g::Grid{T, 3, :cylindrical} = pss.grid
    layout --> (2, 2) 

    cross_section::Symbol, idx::Int = if ismissing(φ) && ismissing(r) && ismissing(z)
        :φ, 1
    elseif !ismissing(φ) && ismissing(r) && ismissing(z)
        φ_rad::T = T(deg2rad(φ))
        while !(g[:φ].interval.left <= φ_rad <= g[:φ].interval.right)
            if φ_rad > g[:φ].interval.right
                φ_rad -= g[:φ].interval.right - g[:φ].interval.left
            elseif φ_rad < g[:φ].interval.left
                φ_rad += g[:φ].interval.right - g[:φ].interval.left
            end
        end
        :φ, searchsortednearest(g[:φ], φ_rad)
    elseif ismissing(φ) && !ismissing(r) && ismissing(z)
        :r, searchsortednearest(g[:r], T(r))
    elseif ismissing(φ) && ismissing(r) && !ismissing(z)
        :z, searchsortednearest(g[:z], T(z))
    else
        error(ArgumentError, ": Only one of the keywords `r, φ, z` is allowed.")
    end
    value::T = if cross_section == :φ
        g[:φ][idx]
    elseif cross_section == :r    
        g[:r][idx]
    elseif cross_section == :z
        g[:z][idx]
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
        ChargeDensity(pss)
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


@recipe function f( pss::PotentialSimulationSetup{T, 3, :cartesian};
                    x = missing,
                    y = 0,
                    z = missing ) where {T}
    g::Grid{T, 3, :cartesian} = pss.grid
    layout --> (2, 2) 

    size --> (1000, 1000)

    cross_section::Symbol, idx::Int = if ismissing(x) && ismissing(y) && ismissing(z)
        :x, searchsortednearest(g[:x], T(0))
    elseif !ismissing(x) && ismissing(y) && ismissing(z)
        :x, searchsortednearest(g[:x], T(x))
    elseif ismissing(x) && !ismissing(y) && ismissing(z)
        :y, searchsortednearest(g[:y], T(y))
    elseif ismissing(x) && ismissing(y) && !ismissing(z)
        :z, searchsortednearest(g[:z], T(z))
    else
        error(ArgumentError, ": Only one of the keywords `x, y, z` is allowed.")
    end

    if cross_section == :x
        x --> g[:x][idx]
    elseif cross_section == :y
        y --> g[:y][idx]
    elseif cross_section == :z
        z --> g[:z][idx]
    end


    @series begin
        subplot := 1
        PointTypes(pss)
    end
    @series begin
        subplot := 2
        ChargeDensity(pss)
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
