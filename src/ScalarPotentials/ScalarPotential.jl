include("EffectiveChargeDensity.jl")
include("DielectricDistribution.jl")
include("ElectricPotential.jl")
include("WeightingPotential.jl")
include("PointTypes.jl")
include("ImpurityScale.jl")

const ScalarPotential{T, N, S, AT} = Union{
    ElectricPotential{T, N, S, AT},
    WeightingPotential{T, N, S, AT},
    PointTypes{T, N, S, AT},
    EffectiveChargeDensity{T, N, S, AT},
    DielectricDistribution{T, N, S, AT},
    ImpurityScale{T, N, S, AT}
}

ScalarPotential(::ElectricPotential, data, grid) = ElectricPotential(data, grid)
ScalarPotential(::WeightingPotential, data, grid) = WeightingPotential(data, grid)
ScalarPotential(::PointTypes, data, grid) = PointTypes(data, grid)
ScalarPotential(::EffectiveChargeDensity, data, grid) = EffectiveChargeDensity(data, grid)
ScalarPotential(::ImpurityScale, data, grid) = ImpurityScale(data, grid)

get_axes_type(::Type{Tuple{AT1, AT2, AT3}}) where {AT1, AT2, AT3} = (AT1, AT2, AT3)

function getindex(ep::P, grid::Grid{T, N, S}) where {T, N, S, P <: ScalarPotential{T, N, S}}
    gridsize::Tuple = size(grid)
    data::Array{T, N} = zeros(T, gridsize)
    ep_itp::Interpolations.Extrapolation{T, N} = interpolated_scalarfield(ep)
    PT = (S == Cylindrical ? CylindricalPoint : CartesianPoint)
    for i1 in eachindex(grid[1])
        for i2 in eachindex(grid[2])
            for i3 in eachindex(grid[3])
                data[i1, i2, i3] = get_interpolation(ep_itp, PT{T}(grid[i1, i2, i3]), S)
            end
        end
    end
    return P(data, grid)
end

function get_2π_potential(sp::ScalarPotential{T, 3, Cylindrical}, axφ::DiscreteAxis{AT, :periodic, :periodic}, int::Interval{:closed, :open, AT}) where {T, AT}
    @assert int.right != 0 "Right boundary of φ interval is not allowed to be 0"
    l::Int = length( axφ )
    Δφ::AT = width(int)
    new_int::Interval{:closed, :open, AT} = Interval{:closed, :open, AT}(0, 2π)
    n::Int = Int(round(T(2π) / Δφ, sigdigits = 6))
    new_ticks::Vector{AT} = Vector{AT}(undef, l * n)
    new_pot = Array{eltype(sp.data), 3}(undef, size(sp, 1), l * n, size(sp, 3))
    for idx_n in 1:n
        for il in 1:l
            idx::Int = il + (idx_n - 1) * l
            new_ticks[idx] = mod((idx_n - 1) * Δφ + axφ.ticks[il], T(2π))
            new_pot[:, idx, :] = sp[:, il, :]
        end
    end
    P = sortperm(new_ticks)
    if new_ticks[P][1] == 0
        new_axφ = DiscreteAxis{AT, :periodic, :periodic}( new_int, new_ticks[P] )
        new_axes = (sp.grid[1], new_axφ, sp.grid[3])
        new_grid = Grid{AT, 3, Cylindrical, typeof(new_axes)}( new_axes )
        ScalarPotential(sp, new_pot[:,P,:], new_grid)
    else
        new_ticks_new_zero::Vector{AT} = Vector{AT}(undef, l * n + 1)
        new_ticks_new_zero[1] = T(0)
        new_ticks_new_zero[2:end] = new_ticks[P]
        new_pot_new_zero = Array{eltype(sp.data), 3}(undef, size(sp, 1), l * n + 1, size(sp, 3))
        new_pot_new_zero[:,1,:] = (eltype(sp.data) <: AbstractFloat) ? sp[:,1,:] + 
                        (sp[:,end,:] - sp[:,1,:]) * (T(0) - new_ticks_new_zero[2]) / (T(2π) - T(new_ticks_new_zero[end]) * n + T(new_ticks_new_zero[2])) :
                        sp[:,1,:]
        new_pot_new_zero[:,2:end,:] = new_pot[:,P,:]
        new_axφ = DiscreteAxis{AT, :periodic, :periodic}( new_int, new_ticks_new_zero )
        new_axes = (sp.grid[1], new_axφ, sp.grid[3])
        new_grid = Grid{AT, 3, Cylindrical, typeof(new_axes)}( new_axes )
        ScalarPotential(sp, new_pot_new_zero, new_grid)
    end    
end



function get_2π_potential(sp::ScalarPotential{T, 3, Cylindrical}, axφ::DiscreteAxis{AT, :reflecting, :reflecting}, int::Interval{:closed, :closed, AT}) where {T, AT}
    @assert int.right != 0 "Right boundary of φ interval is not allowed to be 0"
    l::Int = 2 * length( axφ ) - 2
    Δφ::AT = width(int)
    new_int::Interval{:closed, :open, AT} = Interval{:closed, :open, AT}(int.left, int.right + Δφ)
    new_ticks::Vector{AT} = Vector{AT}(undef, l)
    new_ticks[1:length( axφ )] = collect(axφ.ticks)
    new_pot = Array{eltype(sp.data), 3}(undef, size(sp, 1), l, size(sp, 3))
    new_pot[:, 1:length( axφ ), :] = sp.data[:, :, :]
    for i in 1:(length(axφ) - 2)
        new_ticks[length(axφ) + i] = 2 * int.right - axφ.ticks[length(axφ) - i]
        new_pot[:, length(axφ) + i, :] = sp[:, length(axφ) - i, :]
    end
    new_axφ = DiscreteAxis{AT, :periodic, :periodic}( new_int, new_ticks )
    new_axes = (sp.grid[1], new_axφ, sp.grid[3])
    new_grid = Grid{AT, 3, Cylindrical, typeof(new_axes)}( new_axes )
    ScalarPotential(sp, new_pot, new_grid)
end

function extend_2D_to_3D_by_n_points(sp::ScalarPotential{T, 3, Cylindrical}, axφ::DiscreteAxis{AT, :reflecting, :reflecting}, int::Interval{:closed, :closed, AT},
        n_points_in_φ::Int ) where {T, AT}
    new_int::Interval{:closed, :open, AT} = Interval{:closed, :open, AT}(0, 2π)
    new_ticks::Vector{AT} = Vector{AT}(undef, n_points_in_φ)
    new_pot = Array{eltype(sp.data), 3}(undef, size(sp, 1), n_points_in_φ, size(sp, 3))
    Δφ::AT = 2π / n_points_in_φ
    for i in 1:n_points_in_φ
        new_ticks[i] = (i - 1) * Δφ
        new_pot[:, i, :] = sp[:, 1, :]
    end
    new_axφ = DiscreteAxis{AT, :periodic, :periodic}( new_int, new_ticks )
    new_axes = (sp.grid[1], new_axφ, sp.grid[3])
    new_grid = Grid{AT, 3, Cylindrical, typeof(new_axes)}( new_axes )
    ScalarPotential(sp, new_pot, new_grid)
end

function get_2π_potential(sp::ScalarPotential{T, 3, Cylindrical}; n_points_in_φ::Union{Missing, Int} = missing) where {T}
    axφ::DiscreteAxis{T} = sp.grid[2]
    int::Interval = axφ.interval
    if length(axφ) == 1 # 2D
        if ismissing(n_points_in_φ)
            error(ArgumentError, ": First Argument is a 2D potential (only 1 point in φ). User has to set the keyword `n_points_in_φ::Int` (e.g. `n_points_in_φ = 18`) in order to get a 3D potential.")
        else
            return extend_2D_to_3D_by_n_points(sp, axφ, int, n_points_in_φ)
        end
    else
        if isclosedset(axφ.interval) 
            sp_tmp = get_2π_potential(sp, axφ, int)
            return get_2π_potential(sp_tmp, sp_tmp.grid[2], sp_tmp.grid[2].interval)
        else
            return get_2π_potential(sp, axφ, int)
        end
    end
end