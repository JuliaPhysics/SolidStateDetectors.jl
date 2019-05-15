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


function NamedTuple(ρ::ChargeDensity{T}) where {T <: SSDFloat}
    return (
        grid = NamedTuple(ρ.grid),
        values = ρ.data * internal_voltage_unit,
    )
end
Base.convert(T::Type{NamedTuple}, x::ChargeDensity) = T(x)
    
function ChargeDensity(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(ustrip(nt.values[1]))
    S = get_coordinate_type(grid)
    N = get_number_of_dimensions(grid)
    ChargeDensity{T, N, S}( ustrip.(uconvert.(internal_voltage_unit, nt.values)), grid)
end
Base.convert(T::Type{ChargeDensity}, x::NamedTuple) = T(x)



@recipe function f( ρ::ChargeDensity{T, 3, :cylindrical};
                    r = missing,
                    φ = missing,
                    z = missing ) where {T}
    g::Grid{T, 3, :cylindrical} = ρ.grid
   
    # seriescolor --> :viridis
    st --> :heatmap
    aspect_ratio --> 1
    foreground_color_border --> nothing
    tick_direction --> :out
       
       
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

    @series begin
        if cross_section == :φ
            title --> "Charge Density @$(cross_section) = $(round(rad2deg(value), sigdigits = 2))"
            xlabel --> "r / m"
            ylabel --> "z / m"
            size --> ( 400, 350 / (g[:r][end] - g[:r][1]) * (g[:z][end] - g[:z][1]) )
            g[:r], g[:z], ρ.data[:, idx,:]'
        elseif cross_section == :r
            title --> "Charge Density @$(cross_section) = $(round(value, sigdigits = 2))"
            g[:φ], g[:z], ρ.data[idx,:,:]'
        elseif cross_section == :z
            title --> "Charge Density @$(cross_section) = $(round(value, sigdigits = 2))"
            proj --> :polar
            g[:φ], g[:r], ρ.data[:,:,idx]
        end
    end
end


@recipe function f( ρ::ChargeDensity{T, 3, :cartesian};
                    x = missing,
                    y = missing,
                    z = missing ) where {T}
    g::Grid{T, 3, :cartesian} = ρ.grid
   
    # seriescolor --> :viridis
    st --> :heatmap
    aspect_ratio --> 1
    foreground_color_border --> nothing
    tick_direction --> :out
       
       
    cross_section::Symbol, idx::Int = if ismissing(x) && ismissing(y) && ismissing(z)
        :x, 1
    elseif !ismissing(x) && ismissing(y) && ismissing(z)
        :x, searchsortednearest(g[:x], T(x))
    elseif ismissing(x) && !ismissing(y) && ismissing(z)
        :y, searchsortednearest(g[:y], T(y))
    elseif ismissing(x) && ismissing(y) && !ismissing(z)
        :z, searchsortednearest(g[:z], T(z))
    else
        error(ArgumentError, ": Only one of the keywords `r, y, z` is allowed.")
    end
    value::T = if cross_section == :x
        g[:x][idx]
    elseif cross_section == :y
        g[:y][idx]
    elseif cross_section == :z
        g[:z][idx]
    end


    @series begin
        title --> "Charge Density @$(cross_section) = $(round(value, sigdigits = 2))"
        if cross_section == :x
            xlabel --> "y / m"
            ylabel --> "z / m"
            g[:y], g[:z], ρ.data[idx, :, :]'
        elseif cross_section == :y
            xlabel --> "x / m"
            ylabel --> "z / m"
            g[:x], g[:z], ρ.data[:, idx, :]'
        elseif cross_section == :z
            xlabel --> "x / m"
            ylabel --> "y / m"
            g[:x], g[:y], ρ.data[:,:,idx]'
        end
    end
end

