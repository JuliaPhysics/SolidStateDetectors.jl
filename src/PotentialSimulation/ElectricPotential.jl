struct ElectricPotential{T, N, S} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S}
end

@inline size(ep::ElectricPotential{T, N, S}) where {T, N, S} = size(ep.data)
@inline length(ep::ElectricPotential{T, N, S}) where {T, N, S} = length(ep.data)
@inline getindex(ep::ElectricPotential{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(ep.data, I...)
@inline getindex(ep::ElectricPotential{T, N, S}, i::Int) where {T, N, S} = getindex(ep.data, i)
@inline getindex(ep::ElectricPotential{T, N, S}, s::Symbol) where {T, N, S} = getindex(ep.grid, s)


"""
    ElectricPotential(setup::PotentialSimulationSetup{T, 3, :cylindrical} ; kwargs...)::ElectricPotential{T, 3, :cylindrical}

Extracts the electric potential from `setup` and extrapolate it to an 2π grid.

For 2D grids (r and z) the user has to set the keyword `n_points_in_φ::Int`, e.g.: `n_points_in_φ = 36`.
"""
function ElectricPotential(setup::PotentialSimulationSetup{T, 3, :cylindrical} ; kwargs...)::ElectricPotential{T, 3, :cylindrical} where {T}
    return get_2π_potential(ElectricPotential{T, 3, :cylindrical}(setup.potential, setup.grid); kwargs...)
end

"""
    ElectricPotential(setup::PotentialSimulationSetup{T, 3, :cartesian} ; kwargs...)::ElectricPotential{T, 3, :cartesian}

Extracts the electric potential from `setup`.
"""
function ElectricPotential(setup::PotentialSimulationSetup{T, 3, :cartesian} )::ElectricPotential{T, 3, :cartesian} where {T}
    return ElectricPotential{T, 3, :cartesian}(setup.potential, setup.grid)
end



function NamedTuple(ep::ElectricPotential{T, 3}) where {T}
    return (
        grid = NamedTuple(ep.grid),
        values = ep.data * u"V",
    )
end
Base.convert(T::Type{NamedTuple}, x::ElectricPotential) = T(x)
    
function ElectricPotential(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(ustrip(nt.values[1]))
    S = get_coordinate_system(grid)
    N = get_number_of_dimensions(grid)
    ElectricPotential{T, N, S}( ustrip.(uconvert.(u"V", nt.values)), grid)
end
Base.convert(T::Type{ElectricPotential}, x::NamedTuple) = T(x)


#=
@recipe function f( ep::ElectricPotential{T, 3, :cartesian};
                    # dim = missing, dimvalue = missing,
                    x = missing,
                    y = missing,
                    z = missing,
                    contours_equal_potential=false ) where {T}
    g::Grid{T, 3, :cartesian} = ep.grid
   
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
        title --> "Electric Potential @$(cross_section) = $(round(value, sigdigits = 2))"
        if cross_section == :x
            xlabel --> "y / m"
            ylabel --> "z / m"
            g[:y], g.z, ep.data[idx, :, :]'
        elseif cross_section == :y
            xlabel --> "x / m"
            ylabel --> "z / m"
            g[:x], g.z, ep.data[:, idx, :]'
        elseif cross_section == :z
            xlabel --> "x / m"
            ylabel --> "y / m"
            g[:x], g[:y], ep.data[:,:,idx]'
        end
    end
    if contours_equal_potential
        @series begin
            seriescolor := :thermal
            st := :contours
            if cross_section == :x
                g[:y], g.z, ep.data[idx, :, :]'
            elseif cross_section == :y
                g[:x], g.z, ep.data[:, idx, :]'
            elseif cross_section == :z
                g[:x], g[:y], ep.data[:,:,idx]
            end
        end
    end
end
=#