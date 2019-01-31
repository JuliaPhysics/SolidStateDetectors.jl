
@recipe function f( pss::PotentialSimulationSetup{T, 3, :Cylindrical};
                    r = missing,
                    θ = missing,
                    z = missing,
                    n_points_in_θ = 36 ) where {T}
    g::Grid{T, 3, :Cylindrical} = pss.grid
    layout --> (2, 2) 

    cross_section::Symbol, idx::Int = if ismissing(θ) && ismissing(r) && ismissing(z)
        :θ, 1
    elseif !ismissing(θ) && ismissing(r) && ismissing(z)
        θ_rad::T = T(deg2rad(θ))
        while !(g[:θ].interval.left <= θ_rad <= g[:θ].interval.right)
            if θ_rad > g[:θ].interval.right
                θ_rad -= g[:θ].interval.right - g[:θ].interval.left
            elseif θ_rad < g[:θ].interval.left
                θ_rad += g[:θ].interval.right - g[:θ].interval.left
            end
        end
        :θ, searchsortednearest(g[:θ], θ_rad)
    elseif ismissing(θ) && !ismissing(r) && ismissing(z)
        :r, searchsortednearest(g[:r], T(r))
    elseif ismissing(θ) && ismissing(r) && !ismissing(z)
        :z, searchsortednearest(g[:z], T(z))
    else
        error(ArgumentError, ": Only one of the keywords `r, θ, z` is allowed.")
    end
    value::T = if cross_section == :θ
        g[:θ][idx]
    elseif cross_section == :r    
        g[:r][idx]
    elseif cross_section == :z
        g[:z][idx]
    end

    if cross_section == :θ
        θ --> θ
    elseif cross_section == :z
        z --> z
    elseif cross_section == :r
        r --> r
    end

    @series begin
        subplot := 1
        PointTypes(pss, n_points_in_θ=n_points_in_θ)
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
        ElectricPotential(pss, n_points_in_θ=n_points_in_θ)
    end
end


@recipe function f( pss::PotentialSimulationSetup{T, 3, :Cartesian};
                    x = missing,
                    y = missing,
                    z = missing ) where {T}
    g::Grid{T, 3, :Cartesian} = pss.grid
    layout --> (2, 2) 

    size --> (1000, 1000)

    cross_section::Symbol, idx::Int = if ismissing(x) && ismissing(y) && ismissing(z)
        :x, 1
    elseif !ismissing(x) && ismissing(y) && ismissing(z)
        :x, searchsortednearest(g[:x], T(x))
    elseif ismissing(x) && !ismissing(y) && ismissing(z)
        :y, searchsortednearest(g[:y], T(y))
    elseif ismissing(x) && ismissing(y) && !ismissing(z)
        :z, searchsortednearest(g[:z], T(z))
    else
        error(ArgumentError, ": Only one of the keywords `x, y, z` is allowed.")
    end
    value::T = if cross_section == :x
        g[:x][idx]
    elseif cross_section == :y    
        g[:y][idx]
    elseif cross_section == :z
        g[:z][idx]
    end
    if cross_section == :x
        x --> x
    elseif cross_section == :y
        y --> y
    elseif cross_section == :z
        z --> z
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
