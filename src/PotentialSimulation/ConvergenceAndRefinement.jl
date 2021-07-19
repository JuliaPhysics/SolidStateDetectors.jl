function update_and_get_max_abs_diff!(  fssrb::PotentialSimulationSetupRB{T, N1, N2},
                                        depletion_handling::Val{depletion_handling_enabled}, 
                                        only2d::Val{only_2d} = Val{false}(), 
                                        is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
                                        use_nthreads::Int = Base.Threads.nthreads()
                                        )::T where {T, N1, N2, depletion_handling_enabled, only_2d, _is_weighting_potential}
    tmp_potential::Array{T, N2} = copy(fssrb.potential)
    if depletion_handling_enabled
        update!(fssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
        slopes::Array{T, N2} = tmp_potential - fssrb.potential
        @inbounds for i in 1:19
            tmp_potential[:] = fssrb.potential[:]
            update!(fssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
            slopes += tmp_potential - fssrb.potential
        end
        @inbounds slopes /= 20
        return maximum(abs.(slopes))
    else
        for i in 1:10
            update!(fssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
        end
        max_diff::T = maximum(abs.(tmp_potential - fssrb.potential))
        return max_diff
    end
end

function _update_till_convergence!( fssrb::PotentialSimulationSetupRB{T, N1, N2}, 
                                    convergence_limit::T;
                                    n_iterations_between_checks = 500,
                                    depletion_handling::Val{depletion_handling_enabled} = Val{false}(),
                                    only2d::Val{only_2d} = Val{false}(), 
                                    is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
                                    use_nthreads::Int = Base.Threads.nthreads(), 
                                    max_n_iterations::Int = -1
                                    )::T where {T, N1, N2, depletion_handling_enabled, only_2d, _is_weighting_potential}
    n_iterations::Int = 0
    cf::T = Inf
    cfs::Vector{T} = fill(cf, 10)
    cl::T = _is_weighting_potential ? convergence_limit : abs(convergence_limit * fssrb.bias_voltage) # to get relative change in respect to bias voltage
    # To disable automatically ProgressMeters in CI builds:
    is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")
    prog = ProgressThresh(cl; dt = 0.1, desc = "Convergence: ", output = stderr, enabled = !is_logging(stderr))
    while cf > cl
        for i in 1:n_iterations_between_checks
            update!(fssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
        end
        cf = update_and_get_max_abs_diff!(fssrb, depletion_handling, only2d, is_weighting_potential, use_nthreads)
        @inbounds cfs[1:end-1] = cfs[2:end]
        @inbounds cfs[end] = cf
        slope::T = abs(mean(diff(cfs)))
        ProgressMeter.update!(prog, cf)
        n_iterations += n_iterations_between_checks
        if slope < cl
            # @info "Slope is basically 0 -> Converged: $slope"
            cf = slope
        end
        if max_n_iterations > 0 && n_iterations > max_n_iterations
            # @show n_iterations_between_checks
            @info "Maximum number of iterations reached. (`n_iterations = $(n_iterations)`)"
            break
        end
    end
    if depletion_handling_enabled
        tmp_pointtypes::Array{PointType, N2} = fssrb.pointtypes .& undepleted_bit
        @showprogress "Checking undepleted regions " for i in 1:10
            update!(fssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
            @inbounds for i in eachindex(fssrb.pointtypes)
                if (fssrb.pointtypes[i] & undepleted_bit == 0) && (tmp_pointtypes[i] > 0)
                    fssrb.pointtypes[i] += undepleted_bit
                elseif (fssrb.pointtypes[i] & undepleted_bit > 0) && (tmp_pointtypes[i] == 0)
                    tmp_pointtypes[i] += undepleted_bit
                end
            end
        end
        @inbounds for i in eachindex(fssrb.pointtypes)
            if (fssrb.pointtypes[i] & update_bit == 0)
                fssrb.pointtypes[i] = PointType(0)
            else
                if (fssrb.pointtypes[i] & pn_junction_bit == 0)
                    if fssrb.pointtypes[i] & undepleted_bit > 0 fssrb.pointtypes[i] -= undepleted_bit end
                end
            end
        end
    end

    ProgressMeter.finish!(prog)
    return cf
end


"""
    refine_scalar_potential(p::ScalarPotential{T}, max_diffs::NTuple{3, T}, minimum_distances::NTuple{3, T}; 
        only2d::Val{only_2d} = Val(size(p.data, 2)==1)) where {T, only_2d}

Refine any scalar potential `p`. 

1. Extent the grid to be a closed grid in all dimensions. 
2. Refine the axis of the grid based on `max_diffs` and `minimum_applied_potential`:
   Insert N new ticks between to existing ticks such that the potential difference between each tick becomes
   smaller than `max_diff[i]` (i -> dimension) but that the distances between the ticks stays larger than `minimum_distances[i]`.
3. Create the new data array for the refined grid and fill it by interpolation of the the initial (coarse) grid.
"""
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
    return _convert_closed_potential(typeof(p), ElectricPotential(new_data, new_grid))
end             

function interpolate_closed_potential(p::ScalarPotential, ::Val{true}) where {T}
    interpolate((p.grid.axes[1], p.grid.axes[3]), p.data[:,1,:], Gridded(Linear()))
end
function interpolate_closed_potential(p::ScalarPotential, ::Val{false}) where {T}
    interpolate(p.grid.axes, p.data, Gridded(Linear()))
end

_get_closed_ticks(ticks::Vector{T}, int::ClosedInterval{T}) where {T} = ticks
_get_closed_ticks(ticks::Vector{T}, int::Interval{:closed, :open, T}) where {T} = vcat(ticks, [int.right])
_get_closed_ticks(ticks::Vector{T}, int::Interval{:open, :closed, T}) where {T} = vcat([int.left], ticks)

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

"""
    _get_closed_potential(p::ScalarPotential{T,3,CS}) where {T, CS}

Returnes an closed Grid & Potential:
E.g. if one of the axis is {:closed,:open} it will turn this into {:closed,:closed}
and also extend the `data` field of the potential in the respective dimension and fill
it with the respective values. 
"""
function _get_closed_potential(p::ScalarPotential{T,3,CS}) where {T, CS}
    g = p.grid
    closed_values = _get_closed_values(p.data, 1, g.axes[1].interval)
    closed_values = _get_closed_values(closed_values, 2, g.axes[2].interval)
    closed_values = _get_closed_values(closed_values, 3, g.axes[3].interval)
    closed_axes = broadcast(i -> _get_closed_axis(g.axes[i]), (1, 2, 3))
    AT = typeof(closed_axes)
    closed_grid = Grid{T, 3, CS, AT}(closed_axes)
    ElectricPotential{T, 3, CS, AT}(closed_values, closed_grid) 
end

"""
    _convert_closed_potential(::Type{P}, p::ScalarPotential{T,3,CS}) where {P, T, CS}

Basically the counterpart to `_get_closed_potential`.
"""
function _convert_closed_potential(::Type{P}, p::ScalarPotential{T,3,CS}) where {P, T, CS}
    AT = get_axes_type(P)
    ATs = get_axes_type(AT)
    axs = broadcast(i -> _convert_closed_axis(ATs[i], p.grid.axes[i]), (1, 2, 3))
    values = _reduce_closed_values(ATs, p.data)
    P(values, Grid{T,3,CS,typeof(axs)}(axs))
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
    n_1 = Int.(floor.([maximum(abs.(p.data[i+1,:,:] .- p.data[i,:,:])) for i in 1:size(p.data, 1)-1] ./ max_diffs[1], digits=0)) 
    n_2 = Int.(floor.([maximum(abs.(p.data[:,i+1,:] .- p.data[:,i,:])) for i in 1:size(p.data, 2)-1] ./ max_diffs[2], digits=0)) 
    n_3 = Int.(floor.([maximum(abs.(p.data[:,:,i+1] .- p.data[:,:,i])) for i in 1:size(p.data, 3)-1] ./ max_diffs[3], digits=0)) 
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
    new_axes = broadcast(i -> _refine_axis(p.grid.axes[i], ns[i], sub_widths[i]), (1, 2, 3))
    return typeof(p.grid)(new_axes)
end

function _refine_axis(ax::DiscreteAxis{T, <:Any, <:Any, ClosedInterval{T}}, ns::Vector{Int}, sub_widths::Vector{T}) where {T, I}
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


function _extend_refinement_limits(rl::Union{<:Real, Vector{<:Real}, Tuple{<:Real,<:Real,<:Real}, Vector{<:Tuple{<:Real, <:Real, <:Real}}})
    0
end
_extend_refinement_limits(rl::Real) = (rl, rl, rl )
_extend_refinement_limits(rl::Tuple{<:Real,<:Real,<:Real}) = rl
