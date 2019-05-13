struct WeightingPotential{T, N, S} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S}
end

@inline size(wp::WeightingPotential{T, N, S}) where {T, N, S} = size(wp.data)
@inline length(wp::WeightingPotential{T, N, S}) where {T, N, S} = length(wp.data)
@inline getindex(wp::WeightingPotential{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(wp.data, I...)
@inline getindex(wp::WeightingPotential{T, N, S}, i::Int) where {T, N, S} = getindex(wp.data, i)
@inline getindex(wp::WeightingPotential{T, N, S}, s::Symbol) where {T, N, S} = getindex(wp.grid, s)


function WeightingPotential(fss::PotentialSimulationSetup{T, 3, :cylindrical}; kwargs...)::WeightingPotential{T, 3, :cylindrical} where {T <: SSDFloat}
    return get_2π_potential(WeightingPotential{T, 3, :cylindrical}(fss.potential, fss.grid); kwargs...)
end
function WeightingPotential(fss::PotentialSimulationSetup{T, 3, :cartesian})::WeightingPotential{T, 3, :cartesian} where {T <: SSDFloat}
    return WeightingPotential{T, 3, :cartesian}(fss.potential, fss.grid)
end

@recipe function f( wp::WeightingPotential{T, 3, :cylindrical};
                    r = missing,
                    φ = missing,
                    z = missing ) where {T <: SSDFloat}
    g::Grid{T, 3, :cylindrical} = wp.grid
   
    seriescolor --> :viridis
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
        clims --> (0, 1)
        if cross_section == :φ
            title --> "WeightingPotential Potential @$(cross_section) = $(round(rad2deg(value), sigdigits = 2))"
            xlabel --> "r / m"
            ylabel --> "z / m"
            size --> ( 400, 350 / (g[:r][end] - g[:r][1]) * (g[:z][end] - g[:z][1]) )
            g[:r], g[:z], wp.data[:, idx,:]'
        elseif cross_section == :r
            g[:φ], g[:z], wp.data[idx,:,:]'
        elseif cross_section == :z
            title --> "WeightingPotential Potential @$(cross_section) = $(round(value, sigdigits = 2))"
            proj --> :polar
            g[:φ], g[:r], wp.data[:,:,idx]
        end
    end
end

const ScalarPotential{T, N, S} = Union{ElectricPotential{T, N, S}, WeightingPotential{T, N, S}, PointTypes{T, N, S}, ChargeDensity{T, N, S}}

function get_2π_potential(wp::ScalarPotential{T, 3, :cylindrical}, axφ::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :open, T}) where {T}
    @assert int.right != 0 "Right boundary of φ interval is not allowed to be 0"
    Potential::Type = typeof(wp)
    l::Int = length( axφ )
    Δφ::T = int.right - int.left
    new_int::Interval{:closed, :open, T} = Interval{:closed, :open, T}(0, 2π)
    n::Int = Int(round(T(2π) / Δφ, digits = 5)) 
    new_ticks::Vector{T} = Vector{T}(undef, l * n)
    new_pot::Array{T, 3} = Array{T, 3}(undef, size(wp, 1), l * n, size(wp, 3))
    for idx_n in 1:n
        for il in 1:l
            idx::Int = il + (idx_n - 1) * l
            new_ticks[idx] = (idx_n - 1) * int.right + axφ.ticks[il]
            new_pot[:, idx, :] = wp[:, il, :]
        end
    end
    new_axφ::DiscreteAxis{T, :periodic, :periodic} = DiscreteAxis{T, :periodic, :periodic}( new_int, new_ticks )
    new_grid::Grid{T, 3, :cylindrical} = Grid{T, 3, :cylindrical}( (wp.grid[1], new_axφ, wp.grid[3]) )
    new_wp::Potential = Potential( new_pot, new_grid )
    return new_wp
end

function get_2π_potential(wp::ScalarPotential{T, 3, :cylindrical}, axφ::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval{:closed, :closed, T}) where {T}
    @assert int.right != 0 "Right boundary of φ interval is not allowed to be 0"
    Potential::Type = typeof(wp)
    l::Int = 2 * length( axφ ) - 2
    Δφ::T = int.right - int.left
    new_int::Interval{:closed, :open, T} = Interval{:closed, :open, T}(int.left, int.right + Δφ)
    new_ticks::Vector{T} = Vector{T}(undef, l)
    new_ticks[1:length( axφ )] = collect(axφ.ticks) 
    new_pot::Array{T, 3} = Array{T, 3}(undef, size(wp, 1), l, size(wp, 3))
    new_pot[:, 1:length( axφ ), :] = wp.data[:, :, :]
    for i in 1:(length(axφ) - 2)
        new_ticks[length(axφ) + i] = 2 * int.right - axφ.ticks[length(axφ) - i]
        new_pot[:, length(axφ) + i, :] = wp[:, length(axφ) - i, :]
    end
    new_axφ::DiscreteAxis{T, :periodic, :periodic} = DiscreteAxis{T, :periodic, :periodic}( new_int, new_ticks )
    new_grid::Grid{T, 3, :cylindrical} = Grid{T, 3, :cylindrical}( (wp.grid[1], new_axφ, wp.grid[3]) )
    new_wp::Potential = Potential( new_pot, new_grid )
    return get_2π_potential(new_wp, new_axφ, new_int)
    # return new_wp
end

function extend_2D_to_3D_by_n_points(wp::ScalarPotential{T, 3, :cylindrical}, axφ::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval{:closed, :closed, T},
        n_points_in_φ::Int ) where {T}
    Potential::Type = typeof(wp)
    new_int::Interval{:closed, :open, T} = Interval{:closed, :open, T}(0, 2π)
    new_ticks::Vector{T} = Vector{T}(undef, n_points_in_φ)
    new_pot::Array{T, 3} = Array{T, 3}(undef, size(wp, 1), n_points_in_φ, size(wp, 3))
    Δφ::T = 2π / n_points_in_φ
    for i in 1:n_points_in_φ
        new_ticks[i] = (i - 1) * Δφ
        new_pot[:, i, :] = wp[:, 1, :]
    end
    new_axφ::DiscreteAxis{T, :periodic, :periodic} = DiscreteAxis{T, :periodic, :periodic}( new_int, new_ticks )
    new_grid::Grid{T, 3, :cylindrical} = Grid{T, 3, :cylindrical}( (wp.grid[1], new_axφ, wp.grid[3]) )
    new_wp::Potential = Potential( new_pot, new_grid )
    return new_wp
end

function get_2π_potential(wp::ScalarPotential{T, 3, :cylindrical}; n_points_in_φ::Union{Missing, Int} = missing) where {T}
    axφ::DiscreteAxis{T} = wp.grid[2]
    int::Interval = axφ.interval
    if int.right == 0 && length(axφ) == 1 # 2D
        if ismissing(n_points_in_φ)
            error(ArgumentError, ": First Argument is a 2D potential (only 1 point in φ). User has to set the keyword `n_points_in_φ::Int` (e.g. `n_points_in_φ = 18`) in order to get a 3D potential.")
        else
            return extend_2D_to_3D_by_n_points(wp, axφ, int, n_points_in_φ)
        end
    else
        return get_2π_potential(wp, axφ, int)
    end
end

function WeightingPotential(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(ustrip(nt.values[1]))
    S = get_coordinate_type(grid)
    N = get_number_of_dimensions(grid)
    WeightingPotential{T, N, S}( nt.values, grid)
end
Base.convert(T::Type{WeightingPotential}, x::NamedTuple) = T(x)

function NamedTuple(ep::WeightingPotential{T, 3}) where {T}
    return (
        grid = NamedTuple(ep.grid),
        values = ep.data,
    )
end
Base.convert(T::Type{NamedTuple}, x::WeightingPotential) = T(x)





"""
    PointTypes(setup::PotentialSimulationSetup{T, 3, :cylindrical} ; kwargs...)::PointTypes{T, 3, :cylindrical}

Extracts the electric potential from `setup` and extrapolate it to an 2π grid.

For 2D grids (r and z) the user has to set the keyword `n_points_in_φ::Int`, e.g.: `n_points_in_φ = 36`.
"""
function PointTypes(setup::PotentialSimulationSetup{T, 3, :cylindrical} ; kwargs...)::PointTypes{T, 3, :cylindrical} where {T}
    return get_2π_potential(PointTypes{T, 3, :cylindrical}(setup.pointtypes, setup.grid); kwargs...)
end


@recipe function f( wp::WeightingPotential{T, 3, :cartesian};
                    # dim = missing, dimvalue = missing,
                    x = missing,
                    y = missing,
                    z = missing,
                    contours_equal_potential=false ) where {T}
    g::Grid{T, 3, :cartesian} = wp.grid
   
    seriescolor --> :viridis
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
        title --> "Weighting Potential @$(cross_section) = $(round(value, sigdigits = 2))"
        clims --> (0, 1)
        if cross_section == :x
            xlabel --> "y / m"
            ylabel --> "z / m"
            g[:y], g[:z], wp.data[idx, :, :]'
        elseif cross_section == :y
            xlabel --> "x / m"
            ylabel --> "z / m"
            g[:x], g[:z], wp.data[:, idx, :]'
        elseif cross_section == :z
            xlabel --> "x / m"
            ylabel --> "y / m"
            g[:x], g[:y], wp.data[:,:,idx]'
        end
    end
    if contours_equal_potential
        @series begin
            seriescolor := :thermal
            st := :contours
            if cross_section == :x
                g[:y], g[:z], wp.data[idx, :, :]'
            elseif cross_section == :y
                g[:x], g[:z], wp.data[:, idx, :]'
            elseif cross_section == :z
                g[:x], g[:y], wp.data[:,:,idx]
            end
        end
    end
end
