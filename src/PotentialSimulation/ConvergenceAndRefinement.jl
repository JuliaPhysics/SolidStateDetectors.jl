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
    prog = ProgressThresh(cl, 0.1, "Convergence: ")
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



function _get_refinement_inds( potential::Array{T, 3}, grid::Grid{T, 3, :cylindrical}, 
                max_diffs::NTuple{3, T}, minimum_distances::NTuple{3, T})::NTuple{3, Vector{Int}} where {T <: SSDFloat}
    inds_r::Vector{Int} = Int[]
    inds_φ::Vector{Int} = Int[]
    inds_z::Vector{Int} = Int[]

    ax1::DiscreteAxis{T} = grid[1]
    ax2::DiscreteAxis{T} = grid[2]
    ax3::DiscreteAxis{T} = grid[3]
    int1::Interval = ax1.interval
    int2::Interval = ax2.interval
    int3::Interval = ax3.interval

    refinement_value_r::T = max_diffs[1]
    refinement_value_φ::T = max_diffs[2]
    refinement_value_z::T = max_diffs[3]

    minimum_distance_r::T = minimum_distances[1]
    minimum_distance_φ::T = minimum_distances[2]
    minimum_distance_z::T = minimum_distances[3]

    minimum_distance_φ_deg::T = rad2deg(minimum_distance_φ)

    axr::Vector{T} = collect(grid.axes[1])
    axφ::Vector{T} = collect(grid.axes[2])
    axz::Vector{T} = collect(grid.axes[3])
    cyclic::T = grid.axes[2].interval.right - grid.axes[2].interval.left

    for iz in 1:size(potential, 3)
        z = axz[iz]
        for iφ in 1:size(potential, 2)
            φ = axφ[iφ]
            isum = iz + iφ
            for ir in 1:size(potential, 1)
                r = axr[ir]

                # r
                if ir != size(potential, 1)
                    if abs(potential[ir + 1, iφ, iz] - potential[ir, iφ, iz]) >= refinement_value_r
                        if axr[ir + 1] - axr[ir] > minimum_distance_r
                            if !in(ir, inds_r) push!(inds_r, ir) end
                        end
                    end
                end

                # φ
                if iφ == size(potential, 2)
                    if typeof(int2).parameters[2] == :open # ugly solution for now. should be done over multiple dispatch..
                        if abs(potential[ir, 1, iz] - potential[ir, iφ, iz]) >= refinement_value_φ
                            if cyclic - axφ[iφ] > minimum_distance_φ
                                if !in(iφ, inds_φ) push!(inds_φ, iφ) end
                            end
                        end
                    end
                else
                    if abs(potential[ir, iφ + 1, iz] - potential[ir, iφ, iz]) >= refinement_value_φ
                        if axφ[iφ + 1] - axφ[iφ] > minimum_distance_φ
                            if !in(iφ, inds_φ) push!(inds_φ, iφ) end
                        end
                    end
                end

                # z
                if iz != size(potential, 3)
                    if abs(potential[ir, iφ, iz + 1] - potential[ir, iφ, iz]) >= refinement_value_z
                        if axz[iz + 1] - axz[iz] > minimum_distance_z
                            if !in(iz, inds_z) push!(inds_z, iz) end
                        end
                    end
                end

            end
        end
    end

    if typeof(int2).parameters[2] == :open
        if isodd(length(inds_φ))
            for iφ in 1:size(potential, 2)
                if !in(iφ, inds_φ) 
                    if iφ != size(potential, 2)
                        if axφ[iφ + 1] - axφ[iφ] > minimum_distance_φ
                            append!(inds_φ, iφ)
                            break
                        end
                    else 
                        if cyclic - axφ[iφ] > minimum_distance_φ
                            append!(inds_φ, iφ)
                            break
                        end
                    end
                end
            end
        end
        if isodd(length(inds_φ))
            b::Bool = true
            diffs = diff(axφ)
            while b
                if length(diffs) > 0
                    idx = findmax(diffs)[2]
                    if !in(idx, inds_φ)
                        append!(inds_φ, idx)
                        b = false
                    else
                        deleteat!(diffs, idx)
                    end
                else
                    for i in eachindex(axφ)
                        if !in(i, inds_φ) push!(inds_φ, i) end
                    end
                    b = false
                end
            end
        end    
        @assert iseven(length(inds_φ)) "Refinement would result in uneven grid in φ."
    end

    if isodd(length(inds_z))
        for iz in 1:size(potential, 3)
            if !in(iz, inds_z)
                if iz != size(potential, 3)
                    if axz[iz + 1] - axz[iz] > minimum_distance_z
                        append!(inds_z, iz)
                        break
                    end
                end
            end
        end
    end
    if isodd(length(inds_z)) && ((length(inds_z) + 1) < length(axz))
        b = true
        diffs = diff(axz)
        while b
            if length(diffs) > 0
                idx = findmax(diffs)[2]
                if !in(idx, inds_z)
                    append!(inds_z, idx)
                    b = false
                else
                    deleteat!(diffs, idx)
                end
            else
                for i in eachindex(axz[1:end-1])
                    if !in(i, inds_z) push!(inds_z, i) end
                end
                b = false
            end
        end
    end
    if isodd(length(inds_z)) && (length(inds_z) > 1)
        b = true
        diffs = diff(axz)
        while b && length(diffs) > 1
            idx = findmin(diffs)[2]
            if in(idx, inds_z)
                deleteat!(inds_z, findfirst(i -> i == idx, inds_z))
                b = false
            else
                deleteat!(diffs, idx)
            end
        end
    end
    if isodd(length(inds_z)) inds_z = inds_z[1:end-1] end
    @assert iseven(length(inds_z)) "Refinement would result in uneven grid in z."

    return sort(inds_r), sort(inds_φ), sort(inds_z)
end


function _get_refinement_inds(  potential::Array{T, 3}, grid::Grid{T, 3, :cartesian}, 
                                max_diffs::NTuple{3, T}, minimum_distances::NTuple{3, T})::NTuple{3, Vector{Int}} where {T}
    inds_x::Vector{Int} = Int[]
    inds_y::Vector{Int} = Int[]
    inds_z::Vector{Int} = Int[]

    ax1::DiscreteAxis{T} = grid[1]
    ax2::DiscreteAxis{T} = grid[2]
    ax3::DiscreteAxis{T} = grid[3]
    int1::Interval = ax1.interval
    int2::Interval = ax2.interval
    int3::Interval = ax3.interval

    refinement_value_x::T = max_diffs[1]
    refinement_value_y::T = max_diffs[2]
    refinement_value_z::T = max_diffs[3]

    minimum_distance_x::T = minimum_distances[1]
    minimum_distance_y::T = minimum_distances[2]
    minimum_distance_z::T = minimum_distances[3]

    ax_x::Vector{T} = collect(grid.axes[1])
    ax_y::Vector{T} = collect(grid.axes[2])
    ax_z::Vector{T} = collect(grid.axes[3])

    
    for iz in 1:size(potential, 3)
        z = ax_z[iz]
        for iy in 1:size(potential, 2)
            y = ax_y[iy]
            isum = iz + iy
            for ix in 1:size(potential, 1)
                x = ax_x[ix]

                # x
                if ix != size(potential, 1)
                    if abs(potential[ix + 1, iy, iz] - potential[ix, iy, iz]) >= refinement_value_x
                        if ax_x[ix + 1] - ax_x[ix] > minimum_distance_x
                            if !in(ix, inds_x) push!(inds_x, ix) end
                        end
                    end
                end

                # y
                if iy != size(potential, 2)
                    if abs(potential[ix, iy + 1, iz] - potential[ix, iy, iz]) >= refinement_value_y
                        if ax_y[iy + 1] - ax_y[iy] > minimum_distance_y
                            if !in(iy, inds_y) push!(inds_y, iy) end
                        end
                    end
                end

                # z
                if iz != size(potential, 3)
                    if abs(potential[ix, iy, iz + 1] - potential[ix, iy, iz]) >= refinement_value_z
                        if ax_z[iz + 1] - ax_z[iz] > minimum_distance_z
                            if !in(iz, inds_z) push!(inds_z, iz) end
                        end
                    end
                end

            end
        end
    end

    if isodd(length(inds_x))
        for ix in 1:size(potential, 1)
            if !in(ix, inds_x)
                if ix != size(potential, 1)
                    if ax_x[ix + 1] - ax_x[ix] > minimum_distance_x
                        append!(inds_x, ix)
                        break
                    end
                end
            end
        end
    end
    if isodd(length(inds_x))
        b = true
        diffs = diff(ax_x)
        while b
            if length(diffs) > 0
                idx = findmax(diffs)[2]
                if !in(idx, inds_x)
                    append!(inds_x, idx)
                    b = false
                else
                    deleteat!(diffs, idx)
                end
            else
                for i in eachindex(ax_x[1:end-1])
                    if !in(i, inds_x) push!(inds_x, i) end
                end
                b = false
            end
        end
    end
    if isodd(length(inds_x)) inds_x = inds_x[1:end-1] end
    @assert iseven(length(inds_x)) "Refinement would result in uneven grid in x. This is not allowed since this is the red black dimension."
    

    return sort(inds_x), sort(inds_y), sort(inds_z)
end

function refine(p::ScalarPotential, max_diffs::Tuple{<:Real,<:Real,<:Real}, minimum_distances::Tuple{<:Real,<:Real,<:Real})::typeof(p) 
    T = eltype(p.grid)
    refine_at_inds = _get_refinement_inds(p.data, p.grid, T.(max_diffs), T.(minimum_distances))
    potential, grid = add_points_and_interpolate(p.data, p.grid, refine_at_inds...)
    return typeof(p)(potential, grid)
end