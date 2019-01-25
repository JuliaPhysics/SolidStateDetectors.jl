function update_and_get_max_abs_diff!(  fssrb::PotentialSimulationSetupRB{T, N1, N2},
                                        depletion_handling::Val{depletion_handling_enabled}, only2d::Val{only_2d} = Val{false}(), is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
                                        use_nthreads::Int = Base.Threads.nthreads())::T where {T, N1, N2, depletion_handling_enabled, only_2d, _is_weighting_potential}
    tmp_potential::Array{T, N2} = copy(fssrb.potential)
    update!(fssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
    max_diff::T = maximum(abs.(tmp_potential - fssrb.potential))
    return max_diff
end


function update_till_convergence!(  fssrb::PotentialSimulationSetupRB{T, N1, N2}, 
                                    convergence_limit::Real, bias_voltage::AbstractFloat;
                                    n_iterations_between_checks = 500,
                                    depletion_handling::Val{depletion_handling_enabled} = Val{false}(), only2d::Val{only_2d} = Val{false}(), is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
                                    use_nthreads::Int = Base.Threads.nthreads(), max_n_iterations::Int = -1)::Real where {T, N1, N2, depletion_handling_enabled, only_2d, _is_weighting_potential}
    n_iterations::Int = 0
    cf::T = Inf
    cl::T = abs(convergence_limit * bias_voltage) # to get relative change in respect to bias voltage
    prog = ProgressThresh(cl, 0.1, "Convergence: ")
    while cf > cl
        for i in 1:n_iterations_between_checks
            update!(fssrb, use_nthreads = use_nthreads, depletion_handling = depletion_handling, only2d = only2d, is_weighting_potential = is_weighting_potential)
        end
        cf = update_and_get_max_abs_diff!(fssrb, depletion_handling, only2d, is_weighting_potential, use_nthreads)
        ProgressMeter.update!(prog, cf)
        n_iterations += n_iterations_between_checks
        if max_n_iterations > 0 && n_iterations > max_n_iterations
            @info "Maximum number of iterations reached. (`n_iterations = $(n_iterations)`)"
            break
        end
    end
    ProgressMeter.finish!(prog)
    return cf
end



function check_for_refinement( potential::Array{T, 3}, grid::Grid{T, 3, :Cylindrical}, max_diff::AbstractArray, minimum_distance::AbstractArray)::NTuple{3, Vector{Int}} where {T}
    inds_r::Vector{Int} = Int[]
    inds_θ::Vector{Int} = Int[]
    inds_z::Vector{Int} = Int[]

    ax1::DiscreteAxis{T} = grid[1]
    ax2::DiscreteAxis{T} = grid[2]
    ax3::DiscreteAxis{T} = grid[3]
    int1::Interval = ax1.interval
    int2::Interval = ax2.interval
    int3::Interval = ax3.interval

    refinement_value_r::T = max_diff[1]
    refinement_value_θ::T = max_diff[2]
    refinement_value_z::T = max_diff[3]

    minimum_distance_r::T = minimum_distance[1]
    minimum_distance_θ::T = minimum_distance[2]
    minimum_distance_z::T = minimum_distance[3]

    minimum_distance_θ_deg::T = rad2deg(minimum_distance_θ)

    axr::Vector{T} = collect(grid.axes[1])
    axθ::Vector{T} = collect(grid.axes[2])
    axz::Vector{T} = collect(grid.axes[3])
    cyclic::T = grid.axes[2].interval.right


    for iz in 1:size(potential, 3)
        z = axz[iz]
        for iθ in 1:size(potential, 2)
            θ = axθ[iθ]
            isum = iz + iθ
            for ir in 1:size(potential, 1)
                r = axr[ir]

                # r
                if ir != size(potential, 1)
                    if abs(potential[ir + 1, iθ, iz] - potential[ir, iθ, iz]) >= refinement_value_r
                        if axr[ir + 1] - axr[ir] > minimum_distance_r
                            if !in(ir, inds_r) push!(inds_r, ir) end
                        end
                    end
                end

                # θ
                if iθ == size(potential, 2)
                    if typeof(int2).parameters[2] == :open # ugly solution for now. should be done over multiple dispatch..
                        if abs(potential[ir, 1, iz] - potential[ir, iθ, iz]) >= refinement_value_θ
                            if cyclic - axθ[iθ] > minimum_distance_θ
                                if !in(iθ, inds_θ) push!(inds_θ, iθ) end
                            end
                        end
                    end
                else
                    if abs(potential[ir, iθ + 1, iz] - potential[ir, iθ, iz]) >= refinement_value_θ
                        if axθ[iθ + 1] - axθ[iθ] > minimum_distance_θ
                            if !in(iθ, inds_θ) push!(inds_θ, iθ) end
                        end
                    end
                end

                # z
                if iz != size(potential, 3)
                    if abs(potential[ir, iθ, iz + 1] - potential[ir, iθ, iz]) >= refinement_value_z
                        if axz[iz + 1] - axz[iz] > minimum_distance_z
                            if !in(iz, inds_z) push!(inds_z, iz) end
                        end
                    end
                end

            end
        end
    end

    if typeof(int2).parameters[2] == :open
        if isodd(length(inds_θ))
            for iθ in 1:size(potential, 2)
                if !in(iθ, inds_θ) 
                    if iθ != size(potential, 2)
                        if axθ[iθ + 1] - axθ[iθ] > minimum_distance_θ
                            append!(inds_θ, iθ)
                            break
                        end
                    else 
                        if cyclic - axθ[iθ] > minimum_distance_θ
                            append!(inds_θ, iθ)
                            break
                        end
                    end
                end
            end
        end
        if isodd(length(inds_θ))
            b::Bool = true
            diffs = diff(axθ)
            while b
                if length(diffs) > 0
                    idx = findmax(diffs)[2]
                    if !in(idx, inds_θ)
                        append!(inds_θ, idx)
                        b = false
                    else
                        deleteat!(diffs, idx)
                    end
                else
                    for i in eachindex(axθ)
                        if !in(i, inds_θ) push!(inds_θ, i) end
                    end
                    b = false
                end
            end
        end    
        @assert iseven(length(inds_θ)) "Refinement would result in uneven grid in θ."
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
    if isodd(length(inds_z))
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
    if isodd(length(inds_z)) inds_z = inds_z[1:end-1] end
    @assert iseven(length(inds_z)) "Refinement would result in uneven grid in z."
    

    return sort(inds_r), sort(inds_θ), sort(inds_z)
end
