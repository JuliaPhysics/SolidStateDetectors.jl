struct ChargeDensity{T, N, S} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S}
end

@inline size(ρ::ChargeDensity{T, N, S}) where {T, N, S} = size(ρ.data)
@inline length(ρ::ChargeDensity{T, N, S}) where {T, N, S} = length(ρ.data)
@inline getindex(ρ::ChargeDensity{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(ρ.data, I...)
@inline getindex(ρ::ChargeDensity{T, N, S}, i::Int) where {T, N, S} = getindex(ρ.data, i)
@inline getindex(ρ::ChargeDensity{T, N, S}, s::Symbol) where {T, N, S} = getindex(ρ.grid, s)


function ChargeDensity(fss::PotentialSimulationSetup{T, N, S})::ChargeDensity{T, N, S} where {T, N, S}
    return ChargeDensity{T, N, S}( fss.ρ, fss.grid )
end

@recipe function f( ρ::ChargeDensity{T, 3, :Cylindrical};
                    r = missing,
                    θ = missing,
                    z = missing ) where {T}
    g::Grid{T, 3, :Cylindrical} = ρ.grid
   
    # seriescolor --> :viridis
    st --> :heatmap
    aspect_ratio --> 1
    foreground_color_border --> nothing
    tick_direction --> :out
       
       
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

    @series begin
        if cross_section == :θ
            title --> "Charge Density @$(cross_section) = $(round(rad2deg(value), sigdigits = 2))"
            xlabel --> "r / m"
            ylabel --> "z / m"
            size --> ( 400, 350 / (g[:r][end] - g[:r][1]) * (g[:z][end] - g[:z][1]) )
            g[:r], g[:z], ρ.data[:, idx,:]'
        elseif cross_section == :r
            title --> "Charge Density @$(cross_section) = $(round(value, sigdigits = 2))"
            g[:θ], g[:z], ρ.data[idx,:,:]'
        elseif cross_section == :z
            title --> "Charge Density @$(cross_section) = $(round(value, sigdigits = 2))"
            proj --> :polar
            g[:θ], g[:r], ρ.data[:,:,idx]
        end
    end
end

