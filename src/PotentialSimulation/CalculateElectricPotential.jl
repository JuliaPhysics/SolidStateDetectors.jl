function calculate_electric_potential(  detector::SolidStateDetector{T, :cylindrical};
                                        init_grid_size::NTuple{3, Int} = (10, 10, 10),
                                        init_grid_spacing::Union{Missing, Vector{<:Real}} = missing,
                                        grid::Grid{T, N, S} = Grid(detector, init_grid_size = init_grid_size, init_grid_spacing = init_grid_spacing),
                                        convergence_limit::Real = 5e-6,
                                        max_refinements::Int = 3,
                                        refinement_limits::Vector{<:Real} = [1e-4, 1e-4, 1e-4],
                                        min_grid_spacing::Vector{<:Real} = [1e-4, 1e-2, 1e-4],  # mm, degree, mm
                                        depletion_handling::Bool = false,
                                        use_nthreads::Int = Base.Threads.nthreads(),
                                        sor_consts::Vector{<:Real}=[1.4, 1.85],
                                        max_n_iterations::Int=10000,
                                        verbose::Bool=true,
                                        ) where {T, N, S}
    refinement_limits::Vector{T} = T.(refinement_limits)
    min_grid_spacing::Vector{T} = T.(min_grid_spacing)
    sor_consts::Vector{T} = T.(sor_consts)
    refine::Bool = max_refinements > 0 ? true : false
    only_2d::Bool = length(grid.axes[2]) == 1 ? true : false
    check_grid(grid)
    cyclic::T = grid.axes[2].interval.right - grid.axes[2].interval.left
    n_φ_sym::Int = only_2d ? 1 : Int(round(T(2π) / cyclic, digits = 3))
    if use_nthreads > Base.Threads.nthreads()
        use_nthreads = Base.Threads.nthreads();
        @warn "`use_nthreads` was set to `1`. The environment variable `JULIA_NUM_THREADS` must be set appropriately before the julia session is started."
    end

    n_φ_sym_info_txt = if only_2d
        "φ symmetry: Detector is φ-symmetric -> 2D computation."
    else
        "φ symmetry: calculating just 1/$(n_φ_sym) in φ of the detector."
    end

    fssrb::PotentialSimulationSetupRB{T, 3, 4, :cylindrical} = PotentialSimulationSetupRB(detector, grid);

    if verbose
        println("Electric Potential Calculation\n",
                #"Bulk type: $(detector.bulk_type)\n",
                "Bias voltage: $(fssrb.bias_voltage) V\n",
                "$n_φ_sym_info_txt\n",
                "Precision: $T\n",
                "Convergence limit: $convergence_limit => $(round(abs(fssrb.bias_voltage * convergence_limit), sigdigits=2)) V\n",
                "Threads: $use_nthreads\n",
                "Coordinate system: Cylindrical\n",
                "Initial grid dimension: $(size(grid))\n",
                "Refine? -> $refine\n",
                "Refinement parameters:\n",
                "\tmaximum number of refinements: $(max_refinements)\n",
                "\tminimum grid spacing:\n",
                "\t\tr: $(min_grid_spacing[1]) m\n",
                "\t\tφ: $(min_grid_spacing[2]) rad\n",
                "\t\tz: $(min_grid_spacing[3]) m\n",
                "\tRefinement limits:\n",
                "\t\tr: $(refinement_limits[1]) -> $(round(abs(fssrb.bias_voltage * refinement_limits[1]), sigdigits=2)) V\n",
                "\t\tφ: $(refinement_limits[2]) -> $(round(abs(fssrb.bias_voltage * refinement_limits[2]), sigdigits=2)) V\n",
                "\t\tz: $(refinement_limits[3]) -> $(round(abs(fssrb.bias_voltage * refinement_limits[3]), sigdigits=2)) V\n",
                ""
        )
    end

    n_iterations_between_checks::Int = 500
    update_till_convergence!(   fssrb, convergence_limit, fssrb.bias_voltage,
                                depletion_handling = Val{depletion_handling}(), only2d = Val{only_2d}(),
                                use_nthreads=use_nthreads, max_n_iterations=max_n_iterations, n_iterations_between_checks = n_iterations_between_checks )

    potential::Array{T, 3} = ElectricPotentialArray(fssrb)

    refinement_counter::Int = 0
    if refine
        new_inds = check_for_refinement(potential, grid, abs.(refinement_limits .* fssrb.bias_voltage), min_grid_spacing)
        while !isempty(new_inds[1]) || !isempty(new_inds[2]) || !isempty(new_inds[3])
            refinement_counter += 1
            if refinement_counter <= max_refinements
                potential, grid = add_points_and_interpolate(potential, grid, new_inds...)
                if verbose println("Refinement $refinement_counter:\tNew grid size: $(size(potential))") end
                fssrb = PotentialSimulationSetupRB(detector, grid, potential, sor_consts = refinement_counter == max_refinements ? T[1, 1] : sor_consts)
                n_iterations_between_checks = div(n_iterations_between_checks, 2)
                if n_iterations_between_checks < 50 n_iterations_between_checks = 50 end
                update_till_convergence!(fssrb, convergence_limit, fssrb.bias_voltage, depletion_handling = Val{depletion_handling}(), only2d = Val{only_2d}(), use_nthreads=use_nthreads, max_n_iterations=max_n_iterations, n_iterations_between_checks = n_iterations_between_checks )
                potential = ElectricPotentialArray(fssrb)
                new_inds = check_for_refinement(potential, grid, abs.(refinement_limits .* fssrb.bias_voltage), min_grid_spacing)
            else
                break
            end
        end
    end

    pointtypes::Array{PointType, 3} = PointTypeArray(fssrb)
    # min_or_max_potential_for_bubble_marking::T = fssrb.bulk_is_ptype ? fssrb.minimum_applied_potential + 10 : fssrb.maximum_applied_potential - 10
    @inbounds for i in eachindex(pointtypes)
        if (pointtypes[i] & update_bit == 0)
            pointtypes[i] = PointType(0)
        else
            # if (pointtypes[i] & undepleted_bit > 0) # Deactivated bubble feature for now
            #     if fssrb.bulk_is_ptype
            #         if (min_or_max_potential_for_bubble_marking < potential[i] < fssrb.maximum_applied_potential)
            #             pointtypes[i] += bubble_bit
            #         end
            #     else
            #         if (fssrb.minimum_applied_potential < potential[i] < min_or_max_potential_for_bubble_marking)
            #             pointtypes[i] += bubble_bit
            #         end
            #     end
            # end
            if (pointtypes[i] & pn_junction_bit == 0)
                if pointtypes[i] & undepleted_bit > 0 pointtypes[i] -= undepleted_bit end
                # if pointtypes[i] & bubble_bit > 0 pointtypes[i] -= bubble_bit end
            end
        end
    end

    # Checks
    max::T = maximum(potential)
    min::T = minimum(potential)
    if max > fssrb.maximum_applied_potential
        @warn "Detector not fully depleted at this bias voltage ($(fssrb.bias_voltage) V). At least on grid point has a higher potential value ($(max) V) than the maximum applied potential ($(fssrb.maximum_applied_potential) V). This should not be. However, small overshoots might be due to over relaxation and/or not full convergence."
    end
    if min < fssrb.minimum_applied_potential
        @warn "Detector not fully depleted at this bias voltage ($(fssrb.bias_voltage) V). At least on grid point has a smaller potential value ($(min) V) than the minimum applied potential ($(fssrb.minimum_applied_potential) V). This should not be. However, small overshoots might be due to over relaxation and/or not full convergence."
    end

    return PotentialSimulationSetup{T, N, S}( Grid(fssrb), potential, pointtypes, ChargeDensityArray(fssrb), FixedChargeDensityArray(fssrb), DielektrikumDistributionArray(fssrb)  )
end


"""
    ElectricPotential(detector::SolidStateDetector{T}; kwargs...) where {T}

Calls `calculate_electric_potential(detector; kwargs...)` end extracts only the electric potential and returns it.
"""
function ElectricPotential(detector::SolidStateDetector{T}; kwargs...) where {T}
    setup::PotentialSimulationSetup{T} = calculate_electric_potential(detector; kwargs...)
    return ElectricPotential(setup)
end


function calculate_electric_potential(  detector::SolidStateDetector{T, :cartesian};
                                        init_grid_size::NTuple{3, Int} = (10, 10, 10),
                                        init_grid_spacing::Union{Missing, Vector{<:Real}} = missing,
                                        grid::Grid{T, N, S} = Grid(detector, init_grid_size = init_grid_size, init_grid_spacing = init_grid_spacing),
                                        convergence_limit::Real = 5e-6,
                                        max_refinements::Int = 3,
                                        refinement_limits::Vector{<:Real} = [1e-4, 1e-4, 1e-4],
                                        min_grid_spacing::Vector{<:Real} = [1e-5, 1e-5, 1e-5],  # mm, mm, mm
                                        depletion_handling::Bool = false,
                                        use_nthreads::Int = Base.Threads.nthreads(),
                                        sor_consts::Vector{<:Real}=[1.4],
                                        max_n_iterations::Int=10000,
                                        verbose::Bool=true,
                                        ) where {T, N, S}
    refinement_limits::Vector{T} = T.(refinement_limits)
    min_grid_spacing::Vector{T} = T.(min_grid_spacing)
    sor_consts::Vector{T} = T.(sor_consts)
    refine::Bool = max_refinements > 0 ? true : false

    if use_nthreads > Base.Threads.nthreads()
        use_nthreads = Base.Threads.nthreads();
        @warn "`use_nthreads` was set to `1`. The environment variable `JULIA_NUM_THREADS` must be set appropriately before the julia session is started."
    end

    fssrb::PotentialSimulationSetupRB{T, 3, 4, :cartesian} = PotentialSimulationSetupRB(detector, grid);

    if verbose
        println("Electric Potential Calculation\n",
                #"Bulk type: $(detector.bulk_type)\n",
                "Bias voltage: $(fssrb.bias_voltage) V\n",
                "Precision: $T\n",
                "Convergence limit: $convergence_limit => $(round(abs(fssrb.bias_voltage * convergence_limit), sigdigits=2)) V\n",
                "Threads: $use_nthreads\n",
                "Coordinate system: Cartesian\n",
                "Initial grid dimension: $(size(grid))\n",
                "Refine? -> $refine\n",
                "Refinement parameters:\n",
                "\tmaximum number of refinements: $(max_refinements)\n",
                "\tminimum grid spacing:\n",
                "\t\tx: $(min_grid_spacing[1]) m\n",
                "\t\ty: $(min_grid_spacing[2]) m\n",
                "\t\tz: $(min_grid_spacing[3]) m\n",
                "\tRefinement limits:\n",
                "\t\tx: $(refinement_limits[1]) -> $(round(abs(fssrb.bias_voltage * refinement_limits[1]), sigdigits=2)) V\n",
                "\t\ty: $(refinement_limits[2]) -> $(round(abs(fssrb.bias_voltage * refinement_limits[2]), sigdigits=2)) V\n",
                "\t\tz: $(refinement_limits[3]) -> $(round(abs(fssrb.bias_voltage * refinement_limits[3]), sigdigits=2)) V\n",
                ""
        )
    end

    n_iterations_between_checks::Int = 500
    update_till_convergence!(   fssrb, convergence_limit, fssrb.bias_voltage,
                                depletion_handling = Val{depletion_handling}(),
                                use_nthreads=use_nthreads, max_n_iterations=max_n_iterations, n_iterations_between_checks = n_iterations_between_checks )

    potential::Array{T, 3} = ElectricPotentialArray(fssrb)

    refinement_counter::Int = 0
    if refine
        new_inds = check_for_refinement(potential, grid, abs.(refinement_limits .* fssrb.bias_voltage), min_grid_spacing)
        while !isempty(new_inds[1]) || !isempty(new_inds[2]) || !isempty(new_inds[3])
            refinement_counter += 1
            if refinement_counter <= max_refinements
                potential, grid = add_points_and_interpolate(potential, grid, new_inds...)
                if verbose println("Refinement $refinement_counter:\tNew grid size: $(size(potential))") end
                fssrb = PotentialSimulationSetupRB(detector, grid, potential)
                n_iterations_between_checks = div(n_iterations_between_checks, 2)
                if n_iterations_between_checks < 50 n_iterations_between_checks = 50 end
                update_till_convergence!(fssrb, convergence_limit, fssrb.bias_voltage, depletion_handling = Val{depletion_handling}(), use_nthreads=use_nthreads, max_n_iterations=max_n_iterations, n_iterations_between_checks = n_iterations_between_checks )
                potential = ElectricPotentialArray(fssrb)
                new_inds = check_for_refinement(potential, grid, abs.(refinement_limits .* fssrb.bias_voltage), min_grid_spacing)
            else
                break
            end
        end
    end

    pointtypes::Array{PointType, 3} = PointTypeArray(fssrb)
    # min_or_max_potential_for_bubble_marking::T = fssrb.bulk_is_ptype ? fssrb.minimum_applied_potential + 10 : fssrb.maximum_applied_potential - 10
    @inbounds for i in eachindex(pointtypes)
        if (pointtypes[i] & update_bit == 0)
            pointtypes[i] = PointType(0)
        else
            # if (pointtypes[i] & undepleted_bit > 0) # Deactivated bubble feature for now
            #     if fssrb.bulk_is_ptype
            #         if (min_or_max_potential_for_bubble_marking < potential[i] < fssrb.maximum_applied_potential)
            #             pointtypes[i] += bubble_bit
            #         end
            #     else
            #         if (fssrb.minimum_applied_potential < potential[i] < min_or_max_potential_for_bubble_marking)
            #             pointtypes[i] += bubble_bit
            #         end
            #     end
            # end
            if (pointtypes[i] & pn_junction_bit == 0)
                if pointtypes[i] & undepleted_bit > 0 pointtypes[i] -= undepleted_bit end
                # if pointtypes[i] & bubble_bit > 0 pointtypes[i] -= bubble_bit end
            end
        end
    end

    # Checks
    max::T = maximum(potential)
    min::T = minimum(potential)
    if max > fssrb.maximum_applied_potential
        @warn "Detector not fully depleted at this bias voltage ($(fssrb.bias_voltage) V). At least on grid point has a higher potential value ($(max) V) than the maximum applied potential ($(fssrb.maximum_applied_potential) V). This should not be. However, small overshoots might be due to over relaxation and/or not full convergence."
    end
    if min < fssrb.minimum_applied_potential
        @warn "Detector not fully depleted at this bias voltage ($(fssrb.bias_voltage) V). At least on grid point has a smaller potential value ($(min) V) than the minimum applied potential ($(fssrb.minimum_applied_potential) V). This should not be. However, small overshoots might be due to over relaxation and/or not full convergence."
    end

    return PotentialSimulationSetup{T, N, S}( Grid(fssrb), potential, pointtypes, ChargeDensityArray(fssrb), FixedChargeDensityArray(fssrb), DielektrikumDistributionArray(fssrb)  )
end
