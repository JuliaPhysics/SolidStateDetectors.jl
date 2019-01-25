struct WeightingPotential{T, N, S} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S}
end

@inline size(wp::WeightingPotential{T, N, S}) where {T, N, S} = size(wp.data)
@inline length(wp::WeightingPotential{T, N, S}) where {T, N, S} = length(wp.data)
@inline getindex(wp::WeightingPotential{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(wp.data, I...)
@inline getindex(wp::WeightingPotential{T, N, S}, i::Int) where {T, N, S} = getindex(wp.data, i)
@inline getindex(wp::WeightingPotential{T, N, S}, s::Symbol) where {T, N, S} = getindex(wp.grid, s)


function WeightingPotential(fss::PotentialSimulationSetup{T, N, S}; kwargs...)::WeightingPotential{T, N, S} where {T, N, S}
    return get_2π_potential(WeightingPotential{T, N, S}(fss.potential, fss.grid); kwargs...)
end

@recipe function f( wp::WeightingPotential{T, 3, :Cylindrical};
                    r = missing,
                    θ = missing,
                    z = missing ) where {T}
    g::Grid{T, 3, :Cylindrical} = wp.grid
   
    seriescolor --> :viridis
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
        clims --> (0, 1)
        if cross_section == :θ
            title --> "WeightingPotential Potential @$(cross_section) = $(round(rad2deg(value), sigdigits = 2))"
            xlabel --> "r / m"
            ylabel --> "z / m"
            size --> ( 400, 350 / (g[:r][end] - g[:r][1]) * (g[:z][end] - g[:z][1]) )
            g[:r], g[:z], wp.data[:, idx,:]'
        elseif cross_section == :r
            g[:θ], g[:z], wp.data[idx,:,:]'
        elseif cross_section == :z
            title --> "WeightingPotential Potential @$(cross_section) = $(round(value, sigdigits = 2))"
            proj --> :polar
            g[:θ], g[:r], wp.data[:,:,idx]
        end
    end
end

const ScalarPotential{T, N, S} = Union{ElectricPotential{T, N, S}, WeightingPotential{T, N, S}, PointTypes{T, N, S}, ChargeDensity{T, N, S}}

function get_2π_potential(wp::ScalarPotential{T, 3, :Cylindrical}, axθ::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :open, T}) where {T}
    @assert int.right != 0 "Right boundary of θ interval is not allowed to be 0"
    Potential::Type = typeof(wp)
    l::Int = length( axθ )
    Δθ::T = int.right - int.left
    new_int::Interval{:closed, :open, T} = Interval{:closed, :open, T}(0, 2π)
    n::Int = Int(round(T(2π) / Δθ, digits = 5)) 
    new_ticks::Vector{T} = Vector{T}(undef, l * n)
    new_pot::Array{T, 3} = Array{T, 3}(undef, size(wp, 1), l * n, size(wp, 3))
    for idx_n in 1:n
        for il in 1:l
            idx::Int = il + (idx_n - 1) * l
            new_ticks[idx] = (idx_n - 1) * int.right + axθ.ticks[il]
            new_pot[:, idx, :] = wp[:, il, :]
        end
    end
    new_axθ::DiscreteAxis{T, :periodic, :periodic} = DiscreteAxis{T, :periodic, :periodic}( new_int, new_ticks )
    new_grid::Grid{T, 3, :Cylindrical} = Grid{T, 3, :Cylindrical}( (wp.grid[1], new_axθ, wp.grid[3]) )
    new_wp::Potential = Potential( new_pot, new_grid )
    return new_wp
end

function get_2π_potential(wp::ScalarPotential{T, 3, :Cylindrical}, axθ::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :closed, T}) where {T}
    @assert int.right != 0 "Right boundary of θ interval is not allowed to be 0"
    Potential::Type = typeof(wp)
    l::Int = 2 * length( axθ ) - 2
    Δθ::T = int.right - int.left
    new_int::Interval{:closed, :open, T} = Interval{:closed, :open, T}(int.left, int.right + Δθ)
    new_ticks::Vector{T} = Vector{T}(undef, l)
    new_ticks[1:length( axθ )] = collect(axθ.ticks) 
    new_pot::Array{T, 3} = Array{T, 3}(undef, size(wp, 1), l, size(wp, 3))
    new_pot[:, 1:length( axθ ), :] = wp.data[:, :, :]
    for i in 1:(length(axθ) - 2)
        new_ticks[length(axθ) + i] = 2 * int.right - axθ.ticks[length(axθ) - i]
        new_pot[:, length(axθ) + i, :] = wp[:, length(axθ) - i, :]
    end
    new_axθ::DiscreteAxis{T, :periodic, :periodic} = DiscreteAxis{T, :periodic, :periodic}( new_int, new_ticks )
    new_grid::Grid{T, 3, :Cylindrical} = Grid{T, 3, :Cylindrical}( (wp.grid[1], new_axθ, wp.grid[3]) )
    new_wp::Potential = Potential( new_pot, new_grid )
    return get_2π_potential(new_wp, new_axθ, new_int)
    # return new_wp
end

function extend_2D_to_3D_by_n_points(wp::ScalarPotential{T, 3, :Cylindrical}, axθ::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :closed, T},
        n_points_in_θ::Int ) where {T}
    Potential::Type = typeof(wp)
    new_int::Interval{:closed, :open, T} = Interval{:closed, :open, T}(0, 2π)
    new_ticks::Vector{T} = Vector{T}(undef, n_points_in_θ)
    new_pot::Array{T, 3} = Array{T, 3}(undef, size(wp, 1), n_points_in_θ, size(wp, 3))
    Δθ::T = 2π / n_points_in_θ
    for i in 1:n_points_in_θ
        new_ticks[i] = (i - 1) * Δθ
        new_pot[:, i, :] = wp[:, 1, :]
    end
    new_axθ::DiscreteAxis{T, :periodic, :periodic} = DiscreteAxis{T, :periodic, :periodic}( new_int, new_ticks )
    new_grid::Grid{T, 3, :Cylindrical} = Grid{T, 3, :Cylindrical}( (wp.grid[1], new_axθ, wp.grid[3]) )
    new_wp::Potential = Potential( new_pot, new_grid )
    return new_wp
end

function get_2π_potential(wp::ScalarPotential{T, 3, :Cylindrical}; n_points_in_θ::Union{Missing, Int} = missing) where {T}
    axθ::DiscreteAxis{T} = wp.grid[2]
    int::Interval = axθ.interval
    if int.right == 0 && length(axθ) == 1 # 2D
        if ismissing(n_points_in_θ)
            error(ArgumentError, ": First Argument is a 2D potential (only 1 point in θ). User has to set the keyword `n_points_in_θ::Int` (e.g. `n_points_in_θ = 18`) in order to get a 3D potential.")
        else
            return extend_2D_to_3D_by_n_points(wp, axθ, int, n_points_in_θ)
        end
    else
        return get_2π_potential(wp, axθ, int)
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

function NamedTuple(ep::WeightingPotential{T, 3, :Cylindrical}) where {T}
    return (
        grid = NamedTuple(ep.grid),
        values = ep.data,
    )
end

Base.convert(T::Type{NamedTuple}, x::WeightingPotential) = T(x)

