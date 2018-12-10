function update_and_get_max_abs_diff!(  cgrb::CylindricalRedBlackGrid, pwrb::PrecalculatedWeightsCylindricalRedBlack,
                                        depletion_handling::Val{depletion_handling_enabled}, nthreads::Int=1)::AbstractFloat where {depletion_handling_enabled}
    tmp_potential = copy(cgrb.potential)
    update_pot!(cgrb, pwrb, depletion_handling, nthreads)
    tmp_potential = abs.(tmp_potential - cgrb.potential)
    max_diff = maximum(tmp_potential)
    return max_diff
end


function update_till_convergence!( cgrb::CylindricalRedBlackGrid, pwrb::PrecalculatedWeightsCylindricalRedBlack,
                                   convergence_limit::Real, bias_voltage::AbstractFloat;
                                   n_iterations_between_checks=100,
                                   depletion_handling::Bool=false, nthreads::Int=1, max_n_iterations::Int = -1)::Real
    depletion_handling_enabled = if depletion_handling
        Val{true}()
    else
        Val{false}()
    end

    n_iterations::Int = 0
    cf = update_and_get_max_abs_diff!(cgrb, pwrb, depletion_handling_enabled, nthreads)

    cl = abs(convergence_limit * bias_voltage) # to get relative change in respect to bias voltage
    prog = ProgressThresh(cl, 0.1, "Convergence: ")
    while cf > cl
        for i in 1:n_iterations_between_checks
            update_pot!(cgrb, pwrb, depletion_handling_enabled, nthreads)
        end
        cf = update_and_get_max_abs_diff!(cgrb, pwrb, depletion_handling_enabled, nthreads)
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



function check_for_refinement(  grid::CylindricalGrid, pwrb::PrecalculatedWeightsCylindricalRedBlack, max_diff::AbstractArray, minimum_distance::AbstractArray;
                                info_output::Bool=false)::Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}}
    T = eltype(grid.potential)
    inds_r::Array{Int, 1} = Int[]
    inds_φ::Array{Int, 1} = Int[]
    inds_z::Array{Int, 1} = Int[]

    pt_array = get_PointTypes(grid, pwrb)

    refinement_value_r::T = max_diff[1]
    refinement_value_φ::T = max_diff[2]
    refinement_value_z::T = max_diff[3]

    minimum_distance_r::T = minimum_distance[1]
    minimum_distance_φ::T = minimum_distance[2]
    minimum_distance_z::T = minimum_distance[3]

    minimum_distance_φ_deg::T = deg2rad(minimum_distance_φ)

    if info_output
        println( "Maximum allowed voltage differences:\n\tin r: $(refinement_value_r) V\n\tin φ: $(refinement_value_φ) V\n\tin z: $(refinement_value_z) V")
        println( "Minimum distance between to points: $(T(minimum_distance*1e6)) μm")
        println( "Minimum distance between to points: $(round(minimum_distance_φ_deg, sigdigits=3)) °")
    end


    for iz in 1:size(grid.potential, 3)
        z = grid.z[iz]
        # Δz = z_inv[iz]
        for iφ in 1:size(grid.potential, 2)
            φ = grid.φ[iφ]
            # Δφ = φ_inv[iφ]
            isum = iz + iφ
            for ir in 1:size(grid.potential, 1)
                r = grid.r[ir]
                # Δr = r_inv[ir]

                pt = pt_array[ir, iφ ,iz]

                # r
                if ir != size(grid.potential, 1)
                    if abs(grid.potential[ir + 1, iφ, iz] - grid.potential[ir, iφ, iz]) >= refinement_value_r
                        if grid.r[ir + 1] - grid.r[ir] > minimum_distance_r
                            if !in(ir, inds_r) push!(inds_r, ir) end
                        end
                    end
                end

                # φ
                if iφ == size(grid.potential, 2)
                    if abs(grid.potential[ir, 1, iz] - grid.potential[ir, iφ, iz]) >= refinement_value_φ
                        if grid.cyclic - grid.φ[iφ] > minimum_distance_φ
                            if !in(iφ, inds_φ) push!(inds_φ, iφ) end
                        end
                    end
                else
                    if abs(grid.potential[ir, iφ + 1, iz] - grid.potential[ir, iφ, iz]) >= refinement_value_φ
                        if grid.φ[iφ + 1] - grid.φ[iφ] > minimum_distance_φ
                            if !in(iφ, inds_φ) push!(inds_φ, iφ) end
                        end
                    end
                end

                # z
                if iz != size(grid.potential, 3)
                    if abs(grid.potential[ir, iφ, iz + 1] - grid.potential[ir, iφ, iz]) >= refinement_value_z
                        if grid.z[iz + 1] - grid.z[iz] > minimum_distance_z
                            if !in(iz, inds_z) push!(inds_z, iz) end
                        end
                    end
                end

            end
        end
    end

    if isodd(length(inds_φ))
        for iφ in 1:size(grid.potential, 2)
            if !in(iφ, inds_φ)
                if iφ != size(grid.potential, 2)
                    if grid.φ[iφ + 1] - grid.φ[iφ] > minimum_distance_φ
                        append!(inds_φ, iφ)
                        break
                    end
                else
                    if grid.cyclic - grid.φ[iφ] > minimum_distance_φ
                        append!(inds_φ, iφ)
                        break
                    end
                end
            end
        end
    end
    if isodd(length(inds_φ))
        b::Bool = true
        diffs = diff(grid.φ)
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
                for i in eachindex(grid.φ)
                    if !in(i, inds_φ) push!(inds_φ, i) end
                end
                b = false
            end
        end
    end
    @assert iseven(length(inds_φ)) "Refinement would result in uneven grid in φ."

    if isodd(length(inds_z))
        for iz in 1:size(grid.potential, 3)
            if !in(iz, inds_z)
                if iz != size(grid.potential, 3)
                    if grid.z[iz + 1] - grid.z[iz] > minimum_distance_z
                        append!(inds_z, iz)
                        break
                    end
                end
            end
        end
    end
    if isodd(length(inds_z))
        b = true
        diffs = diff(grid.z)
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
                for i in eachindex(grid.z[1:end-1])
                    if !in(i, inds_z) push!(inds_z, i) end
                end
                b = false
            end
        end
    end
    if isodd(length(inds_z)) inds_z = inds_z[1:end-1] end
    @assert iseven(length(inds_z)) "Refinement would result in uneven grid in z."

    return sort(inds_r), sort(inds_φ), sort(inds_z)
end
