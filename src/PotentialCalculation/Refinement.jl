# """
#     refine_scalar_potential(p::ScalarPotential{T}, max_diffs::NTuple{3, T}, minimum_distances::NTuple{3, T}; 
#         only2d::Val{only_2d} = Val(size(p.data, 2)==1)) where {T, only_2d}
# 
# Refine any scalar potential `p`. 
# 
# 1. Extent the grid to be a closed grid in all dimensions. 
# 2. Refine the axis of the grid based on `max_diffs` and `minimum_applied_potential`:
#    Insert N new ticks between to existing ticks such that the potential difference between each tick becomes
#    smaller than `max_diff[i]` (i -> dimension) but that the distances between the ticks stays larger than `minimum_distances[i]`.
# 3. Create the new data array for the refined grid and fill it by interpolation of the the initial (coarse) grid.
# """
function refine_scalar_potential(p::ScalarPotential{T}, max_diffs::NTuple{3, T}, minimum_distances::NTuple{3, T}; 
        only2d::Val{only_2d} = Val(size(p.data, 2)==1)) where {T, only_2d}
    closed_potential = _get_closed_potential(p)
    new_grid = _create_refined_grid(closed_potential, T.(max_diffs), T.(minimum_distances))
    new_data = Array{T, 3}(undef, size(new_grid))
    if only_2d
        int = interpolate_closed_potential(closed_potential, only2d)
        for i3 in axes(new_data, 3)
            x3 = new_grid.axes[3].ticks[i3]
            for i1 in axes(new_data, 1)
                x1 = new_grid.axes[1].ticks[i1]
                new_data[i1, 1, i3] = int(x1, x3)
            end
        end
    else
        int = interpolate_closed_potential(closed_potential, only2d)
        for i3 in axes(new_data, 3)
            x3 = new_grid.axes[3].ticks[i3]
            for i2 in axes(new_data, 2)
                x2 = new_grid.axes[2].ticks[i2]
                for i1 in axes(new_data, 1)
                    x1 = new_grid.axes[1].ticks[i1]
                    new_data[i1, i2, i3] = int(x1, x2, x3)
                end
            end
        end
    end
    return _convert_to_original_potential(p, new_data, new_grid)
end             

function interpolate_closed_potential(p::ScalarPotential{T}, ::Val{true}) where {T}
    interpolate!((p.grid.axes[1], p.grid.axes[3]), p.data[:,1,:], Gridded(Linear()))
end
function interpolate_closed_potential(p::ScalarPotential{T}, ::Val{false}) where {T}
    interpolate!(p.grid.axes, p.data, Gridded(Linear()))
end

_get_closed_ticks(ticks::Vector{T}, int::ClosedInterval{T}) where {T} = ticks
_get_closed_ticks(ticks::Vector{T}, int::Interval{:closed, :open, T}) where {T} = vcat(ticks, int.right)
_get_closed_ticks(ticks::Vector{T}, int::Interval{:open, :closed, T}) where {T} = vcat(int.left, ticks)

_get_closed_values(values::Array{T,3}, dim::Int, int::ClosedInterval{T}) where {T} = values

function _get_closed_values(values::Array{T,3}, dim::Int, int::Interval{:closed, :open, T})::Array{T,3} where {T} 
    if dim == 1
        cat(values, values[1:1,:,:], dims=dim)
    elseif dim == 2
        cat(values, values[:,1:1,:], dims=dim)
    else # dim == 3
        cat(values, values[:,:,1:1], dims=dim)
    end
end
function _get_closed_values(values::Array{T,3}, dim::Int, int::Interval{:open, :closed, T})::Array{T,3} where {T} 
    if dim == 1
        cat(values[end:end,:,:], values, dims=dim)
    elseif dim == 2
        cat(values[:,end:end,:], values, dims=dim)
    else # dim == 3
        cat(values[:,:,end:end], values, dims=dim)
    end
end

function _get_closed_axis(ax::DiscreteAxis{T, BL, BR}) where {T, BL, BR}
    ticks = _get_closed_ticks(ax.ticks, ax.interval)
    return DiscreteAxis{T, BL, BR, ClosedInterval{T}}(ClosedInterval(ax.interval.left, ax.interval.right), ticks)
end

# """
#     _get_closed_potential(p::ScalarPotential{T,3,CS}) where {T, CS}
# 
# Returns an closed Grid & Potential:
# E.g. if one of the axis is {:closed,:open} it will turn this into {:closed,:closed}
# and also extend the `data` field of the potential in the respective dimension and fill
# it with the respective values.
# """
function _get_closed_potential(p::ScalarPotential{T,3,CS}) where {T, CS}
    grid = p.grid
    closed_values = _get_closed_values(p.data, 1, grid.axes[1].interval)
    closed_values = _get_closed_values(closed_values, 2, grid.axes[2].interval)
    closed_values = _get_closed_values(closed_values, 3, grid.axes[3].interval)
    closed_axes = broadcast(i -> _get_closed_axis(grid.axes[i]), (1, 2, 3))
    AT = typeof(closed_axes)
    closed_grid = Grid{T, 3, CS, AT}(closed_axes)
    ScalarPotential(p, closed_values, closed_grid)
end

# """
#     _convert_to_original_potential(::Type{P}, data, grid) where {P, T, CS}
# 
# Basically the counterpart to `_get_closed_potential`.
# """
function _convert_to_original_potential(p::ScalarPotential{T,3,CS,ATO}, data::Array{T, 3}, grid::Grid{T, 3, CS, AT}) where {T, CS, AT, ATO}
    ATs = get_axes_type(ATO)
    axs = broadcast(i -> _convert_closed_axis(ATs[i], grid.axes[i]), (1, 2, 3))
    new_data = _reduce_closed_values(ATs, data)
    new_grid = Grid{T,3,CS,ATO}(axs)
    ScalarPotential(p, new_data, new_grid)
end

_convert_closed_axis(::Type{DiscreteAxis{T, BL, BR, I}}, ax::DiscreteAxis{T, BL, BR, I}) where {T, BL, BR, I} = ax

function _convert_closed_axis(::Type{DiscreteAxis{T, BL, BR, Interval{:closed, :open, T}}}, ax::DiscreteAxis{T, BL, BR, ClosedInterval{T}}) where {T, BL, BR}
    DiscreteAxis{T, BL, BR, Interval{:closed, :open, T}}(Interval{:closed, :open, T}(ax.interval.left, ax.interval.right), ax.ticks[1:end-1])
end
function _convert_closed_axis(::Type{DiscreteAxis{T, BL, BR, Interval{:open, :closed, T}}}, ax::DiscreteAxis{T, BL, BR, ClosedInterval{T}}) where {T, BL, BR}
    DiscreteAxis{T, BL, BR, Interval{:open, :closed, T}}(Interval{:open, :closed, T}(ax.interval.left, ax.interval.right), ax.ticks[2:end])
end

function _reduce_closed_values(axTypes::Tuple, a::Array{T, 3}) where {T} 
   inds = broadcast(i -> _get_ind_range_of_closed_axis(size(a, i), axTypes[i]), (1, 2, 3))
   a[inds...]
end

_get_ind_range_of_closed_axis(l::Int, ::Type{DiscreteAxis{T, BL, BR, ClosedInterval{T}}}) where {T, BL, BR} = 1:l
_get_ind_range_of_closed_axis(l::Int, ::Type{DiscreteAxis{T, BL, BR, Interval{:closed, :open, T}}}) where {T, BL, BR} = 1:(l-1)
_get_ind_range_of_closed_axis(l::Int, ::Type{DiscreteAxis{T, BL, BR, Interval{:open, :closed, T}}}) where {T, BL, BR} = 2:l


function _create_refined_grid(p::ScalarPotential{T,3}, max_diffs::NTuple{3, T}, minimum_distances::NTuple{3, T}) where {T}
    max_diffs = broadcast(md -> iszero(md) ? Inf : md, max_diffs)
    n_1 = floor.(Int, [maximum(abs.(p.data[i+1,:,:] .- p.data[i,:,:])) for i in 1:size(p.data, 1)-1] ./ max_diffs[1]) 
    n_2 = floor.(Int, [maximum(abs.(p.data[:,i+1,:] .- p.data[:,i,:])) for i in 1:size(p.data, 2)-1] ./ max_diffs[2]) 
    n_3 = floor.(Int, [maximum(abs.(p.data[:,:,i+1] .- p.data[:,:,i])) for i in 1:size(p.data, 3)-1] ./ max_diffs[3]) 
    ns = (n_1, n_2, n_3)
    widths = diff.((p.grid.axes[1].ticks, p.grid.axes[2].ticks, p.grid.axes[3].ticks))
    sub_widths = broadcast(ia -> [widths[ia][i] / (ns[ia][i]+1) for i in eachindex(ns[ia])], (1,2,3))
    for ia in 1:3
        for i in eachindex(ns[ia])
            while sub_widths[ia][i] < minimum_distances[ia] && ns[ia][i] > 0
                ns[ia][i] -= 1
                sub_widths[ia][i] = widths[ia][i] / (ns[ia][i]+1)
            end
        end
    end
    for i in 1:3 # always add an even number of ticks
        if isodd(sum(ns[i])) 
            i_max_width = findmax(sub_widths[i])[2]
            ns[i][i_max_width] += 1 
            sub_widths[i][i_max_width] = widths[i][i_max_width] / (ns[i][i_max_width]+1)
        end
    end
    new_axes = broadcast(i -> _refine_axis(p.grid.axes[i], ns[i], sub_widths[i]), (1, 2, 3))
    return typeof(p.grid)(new_axes)
end

function _refine_axis(ax::DiscreteAxis{T, <:Any, <:Any, ClosedInterval{T}}, ns::Vector{Int}, sub_widths::Vector{T}) where {T}
    @assert length(ns) == length(ax.ticks)-1 # for ClosedInterval axis
    ticks = Vector{T}(undef, length(ax.ticks) + sum(ns))
    i = 1 
    for j in eachindex(ns)
        ticks[i] = ax.ticks[j]
        for k in 1:ns[j]
            i += 1 
            ticks[i] = ax.ticks[j] + k * sub_widths[j]
        end
        i += 1 
    end
    ticks[end] = ax.ticks[end]
    typeof(ax)(ax.interval, ticks)
end


_extend_refinement_limits(rl::Real) = (rl, rl, rl )
_extend_refinement_limits(rl::Tuple{<:Real,<:Real,<:Real}) = rl

@inline function has_surface_points(slice::AbstractArray{PointType})::Bool
    return any(is_in_inactive_layer, slice)
end

function _refine_axis_surface( ax::DiscreteAxis{T, <:Any, <:Any, ClosedInterval{T}}, surface_intervals::AbstractVector{Bool}, min_spacing::T;
    extra_before::Int = 5,  # intervals to refine before first surface interval
    extra_after::Int = 5    # intervals to refine after last surface interval
) where {T}

    old_ticks = ax.ticks
    n_int = length(surface_intervals)

    # Find first and last surface intervals
    first_surface = findfirst(surface_intervals)
    last_surface  = findlast(surface_intervals)

    if first_surface === nothing || last_surface === nothing
        return ax # Nothing to refine
    end

    # Create flags for all intervals to refine
    refine_flags = copy(surface_intervals)

    # Add extra intervals before first surface interval
    start_idx = max(first_surface - extra_before, 1)
    refine_flags[start_idx:first_surface-1] .= true

    # Add extra intervals after last surface interval
    end_idx = min(last_surface + extra_after, n_int)
    refine_flags[last_surface+1:end_idx] .= true

    # Merge consecutive intervals to refine
    merged = Vector{UnitRange{Int}}()
    i = 1
    while i <= n_int
        if refine_flags[i]
            j = i
            while j < n_int && refine_flags[j+1]
                j += 1
            end
            push!(merged, i:j)
            i = j + 1
        else
            i += 1
        end
    end

    # Compute number of points to add per interval
    ns = zeros(Int, n_int)
    for r in merged
        for i in r
            Δ = old_ticks[i+1] - old_ticks[i]
            if Δ > min_spacing
                ns[i] = ceil(Int, Δ/min_spacing) - 1
            end
        end
    end

    # Subdivide intervals
    sub_widths = [(old_ticks[i+1] - old_ticks[i]) / (ns[i]+1) for i in 1:n_int]
    ticks = Vector{T}(undef, length(old_ticks) + sum(ns))
    i_tick = 1
    for j in 1:n_int
        ticks[i_tick] = old_ticks[j]
        for k in 1:ns[j]
            i_tick += 1
            ticks[i_tick] = old_ticks[j] + k*sub_widths[j]
        end
        i_tick += 1
    end
    ticks[end] = old_ticks[end]

    return typeof(ax)(ax.interval, ticks)
end
