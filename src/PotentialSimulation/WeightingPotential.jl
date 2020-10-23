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

function WeightingPotential(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(ustrip(nt.values[1]))
    S = get_coordinate_system(grid)
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

#=
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
        title --> "Weighting Potential @$(cross_section) = $(round(value, sigdigits = 2))"
        clims --> (0, 1)
        if cross_section == :x
            xlabel --> "y / m"
            ylabel --> "z / m"
            g[:y], g.z, wp.data[idx, :, :]'
        elseif cross_section == :y
            xlabel --> "x / m"
            ylabel --> "z / m"
            g[:x], g.z, wp.data[:, idx, :]'
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
                g[:y], g.z, wp.data[idx, :, :]'
            elseif cross_section == :y
                g[:x], g.z, wp.data[:, idx, :]'
            elseif cross_section == :z
                g[:x], g[:y], wp.data[:,:,idx]
            end
        end
    end
end
=#