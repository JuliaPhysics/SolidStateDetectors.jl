const ScalarPotential{T, N, S, AT} = Union{ElectricPotential{T, N, S, AT}, WeightingPotential{T, N, S, AT}, PointTypes{T, N, S, AT}, EffectiveChargeDensity{T, N, S, AT}}

get_axes_type(p::ScalarPotential{T, N, S, AT}) where {T, N, S, AT} = AT
get_axes_type(::Type{<:ScalarPotential{T, N, S, AT}}) where {T, N, S, AT} = AT
get_axes_type(::Type{Tuple{AT1, AT2, AT3}}) where {AT1, AT2, AT3} = (AT1, AT2, AT3)

function getindex(ep::P, g::Grid{T, N, S}) where {T, N, S, P <: ScalarPotential{T, N, S}}
    gridsize::Tuple = size(g)
    data::Array{T, N} = zeros(T, gridsize)
    ep_itp::Interpolations.Extrapolation{T, N} = interpolated_scalarfield(ep)
    point = (S == Cylindrical ? CylindricalPoint : CartesianPoint)
    for i1 in eachindex(g[1])
        for i2 in eachindex(g[2])
            for i3 in eachindex(g[3])
                data[i1, i2, i3] = get_interpolation(ep_itp, point{T}(g[i1, i2, i3]), S)
            end
        end
    end
    return P(data, g)
end


function get_2π_potential(wp::ScalarPotential{T, 3, Cylindrical}, axφ::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :open, T}) where {T}
    @assert int.right != 0 "Right boundary of φ interval is not allowed to be 0"
    Potential::Type = typeof(wp)
    l::Int = length( axφ )
    Δφ::T = int.right - int.left
    new_int::Interval{:closed, :open, T} = Interval{:closed, :open, T}(0, 2π)
    n::Int = Int(round(T(2π) / Δφ, sigdigits = 6))
    new_ticks::Vector{T} = Vector{T}(undef, l * n)
    new_pot::Array{T, 3} = Array{T, 3}(undef, size(wp, 1), l * n, size(wp, 3))
    for idx_n in 1:n
        for il in 1:l
            idx::Int = il + (idx_n - 1) * l
            new_ticks[idx] = (idx_n - 1) * Δφ + axφ.ticks[il]
            new_pot[:, idx, :] = wp[:, il, :]
        end
    end
    new_axφ::DiscreteAxis{T, :periodic, :periodic} = DiscreteAxis{T, :periodic, :periodic}( new_int, new_ticks )
    new_grid::Grid{T, 3, Cylindrical} = Grid{T, 3, Cylindrical}( (wp.grid[1], new_axφ, wp.grid[3]) )
    new_wp::Potential = Potential( new_pot, new_grid )
    return new_wp
end

function get_2π_potential(wp::ScalarPotential{T, 3, Cylindrical}, axφ::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval{:closed, :closed, T}) where {T}
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
    new_grid::Grid{T, 3, Cylindrical} = Grid{T, 3, Cylindrical}( (wp.grid[1], new_axφ, wp.grid[3]) )
    new_wp::Potential = Potential( new_pot, new_grid )
    return get_2π_potential(new_wp, new_axφ, new_int)
    # return new_wp
end

function extend_2D_to_3D_by_n_points(wp::ScalarPotential{T, 3, Cylindrical}, axφ::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval{:closed, :closed, T},
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
    new_grid::Grid{T, 3, Cylindrical} = Grid{T, 3, Cylindrical}( (wp.grid[1], new_axφ, wp.grid[3]) )
    new_wp::Potential = Potential( new_pot, new_grid )
    return new_wp
end

function get_2π_potential(wp::ScalarPotential{T, 3, Cylindrical}; n_points_in_φ::Union{Missing, Int} = missing) where {T}
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