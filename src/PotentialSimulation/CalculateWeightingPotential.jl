"""
    calculate_weighting_potential(det::SolidStateDetector{T}, channel_id::Int; <keyword arguments>) where {T}


Compute the weighting potential for the given Detector `det` on an adaptive grid
through successive over relaxation. It returns a collection struct `PotentialSimulationSetup{T}` which stores
the potential, the charge density, the dielectric distribution, pointtypes and the final grid.

There are serveral `<keyword arguments>` which can be used to tune the computation:

# Keywords
- `convergence_limit::Real`: `convergence_limit` times the bias voltage sets the convergence limit of the relaxation. The convergence value is the absolute maximum difference of the potential between two iterations of all grid points. Default of `convergence_limit` is `5e-6` (times bias voltage).
- `max_refinements::Int`: number of maximum refinements. Default is `2`. Set it to `0` to switch off refinement.
- `refinement_limits::Vector{Real}`: vector of refinement limits for each dimension (in case of cylindrical coordinates the order is `r`, `φ`, `z`). A refinement limit (e.g. `refinement_limits[1]`) times the bias voltage of the detector `det` is the maximum allowed voltage difference between two neighbouring grid points in the respective dimension. When the difference is larger, new points are created inbetween. Default is `[1e-4, 1e-4, 1e-4]`.
- `init_grid_spacing::Vector{Real}`: vector of the initial distances between two grid points for each dimension. For normal coordinates the unit is meter. For angular coordinates, the unit is radiance. It prevents the refinement to make the grid to fine. Default is `[0.005, 10.0, 0.005]``.
- `min_grid_spacing::Vector{Real}`: vector of the mimimum allowed distance between two grid points for each dimension. For normal coordinates the unit is meter. For angular coordinates, the unit is radiance. It prevents the refinement to make the grid to fine. Default is [`1e-4`, `1e-2`, `1e-4`].
- `grid::Grid{T, N, S}`: Initial grid used to start the simulation. Default is `Grid(detector, init_grid_spacing=init_grid_spacing)`.
- `depletion_handling::Bool`: enables the handling of undepleted regions. Default is false.
- `use_nthreads::Int`: Number of threads to use in the computation. Default is `Base.Threads.nthreads()`. The environment variable `JULIA_NUM_THREADS` must be set appropriately before the Julia session was started (e.g. `export JULIA_NUM_THREADS=8` in case of bash).
- `sor_consts::Vector{<:Real}`: Two element array. First element contains the SOR constant for `r` = 0. Second contains the constant at the outer most grid point in `r`. A linear scaling is applied in between. First element should be smaller than the second one and both should be ∈ [1.0, 2.0]. Default is [1.4, 1.85].
- `max_n_iterations::Int`: Set the maximum number of iterations which are performed after each grid refinement. Default is `10000`. If set to `-1` there will be no limit.
- `verbose::Bool=true`: Boolean whether info output is produced or not.
- `init_grid_spacing::Vector{<:Real}`: Initial spacing of the grid. Default is [2e-3, 5, 2e-3] <=> [2mm, 5 degree, 2mm ]
"""
function calculate_weighting_potential( detector::SolidStateDetector{T, :cylindrical}, channel_id::Int;
                                        init_grid_size::NTuple{3, Int} = (10, 10, 10),
                                        init_grid_spacing::Union{Missing, Vector{<:Real}} = missing,
                                        grid::Grid{T, N, S} = Grid(detector, init_grid_size = init_grid_size, init_grid_spacing = init_grid_spacing,
                                                                    for_weighting_potential = true),
                                        convergence_limit::Real = 5e-6,
                                        max_refinements::Int = 3,
                                        refinement_limits::Vector{<:Real} = [1e-4, 1e-3, 1e-4],
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

    fssrb::PotentialSimulationSetupRB{T, 3, 4, :cylindrical} = PotentialSimulationSetupRB(detector, grid, weighting_potential_contact_id = channel_id);

    if verbose
        println("Weighting Potential Calculation\n",
                "$n_φ_sym_info_txt\n",
                "Precision: $T\n",
                "Convergence limit: $convergence_limit => $(round(abs(convergence_limit), sigdigits=2)) V\n",
                "Threads: $use_nthreads\n",
                "Initial grid dimension: $(size(grid))\n",
                "Refine? -> $refine\n",
                "Refinement parameters:\n",
                "\tmaximum number of refinements: $(max_refinements)\n",
                "\tminimum grid spacing:\n",
                "\t\tr: $(min_grid_spacing[1]) m\n",
                "\t\tφ: $(min_grid_spacing[2]) rad\n",
                "\t\tz: $(min_grid_spacing[3]) m\n",
                "\tRefinement limits:\n",
                "\t\tr: $(refinement_limits[1]) -> $(round(abs(refinement_limits[1]), sigdigits=2)) V\n",
                "\t\tφ: $(refinement_limits[2]) -> $(round(abs(refinement_limits[2]), sigdigits=2)) V\n",
                "\t\tz: $(refinement_limits[3]) -> $(round(abs(refinement_limits[3]), sigdigits=2)) V\n",
                ""
        )
    end

    update_till_convergence!(   fssrb, convergence_limit, T(1),
                                depletion_handling = Val{false}(), only2d = Val{only_2d}(), is_weighting_potential = Val{true}(),
                                use_nthreads=use_nthreads, max_n_iterations=max_n_iterations    )

    potential::Array{T, 3} = ElectricPotentialArray(fssrb)

    refinement_counter::Int = 0
    if refine
        new_inds = check_for_refinement(potential, grid, abs.(refinement_limits), min_grid_spacing)
        while !isempty(new_inds[1]) || !isempty(new_inds[2]) || !isempty(new_inds[3])
            refinement_counter += 1
            if refinement_counter <= max_refinements
                potential, grid = add_points_and_interpolate(potential, grid, new_inds...)
                if verbose println("Refinement $refinement_counter:\tNew grid size: $(size(potential))") end
                fssrb = PotentialSimulationSetupRB(detector, grid, potential, weighting_potential_contact_id = channel_id)
                update_till_convergence!(   fssrb, convergence_limit, T(1),
                depletion_handling = Val{false}(), only2d = Val{only_2d}(), is_weighting_potential = Val{true}(),
                use_nthreads=use_nthreads, max_n_iterations=max_n_iterations)
                potential = ElectricPotentialArray(fssrb)
                new_inds = check_for_refinement(potential, grid, abs.(refinement_limits), min_grid_spacing)
            else
                break
            end
        end
    end

    return PotentialSimulationSetup{T, N, S}( Grid(fssrb), potential, PointTypeArray(fssrb), ChargeDensityArray(fssrb), DielektrikumDistributionArray(fssrb)  )
end

"""
    WeightingPotential(detector::SolidStateDetector{T}, channel_id::Int; kwargs...) where {T}

Calls `calculate_weighting_potential(detector, channel_id; kwargs...)` end extracts only the weighting potential and returns it.
"""
function WeightingPotential(detector::SolidStateDetector{T}, channel_id::Int; kwargs...) where {T}
    setup::PotentialSimulationSetup{T} = calculate_weighting_potential(detector, channel_id; kwargs...)
    return WeightingPotential(setup)
end



function calculate_weighting_potential( detector::SolidStateDetector{T, :cartesian}, channel_id::Int;
                                        init_grid_size::NTuple{3, Int} = (10, 10, 10),
                                        init_grid_spacing::Union{Missing, Vector{<:Real}} = missing,
                                        grid::Grid{T, N, S} = Grid(detector, init_grid_size = init_grid_size, init_grid_spacing = init_grid_spacing),
                                        convergence_limit::Real = 5e-6,
                                        max_refinements::Int = 3,
                                        refinement_limits::Vector{<:Real} = [1e-5, 1e-5, 1e-5],
                                        min_grid_spacing::Vector{<:Real} = [1e-6, 1e-6, 1e-6],  
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

    fssrb::PotentialSimulationSetupRB{T, 3, 4, :cartesian} = PotentialSimulationSetupRB(detector, grid, weighting_potential_contact_id = channel_id);

    if verbose
        println("WeightingPotential Potential Calculation\n",
                "Precision: $T\n",
                "Convergence limit: $convergence_limit => $(round(abs(convergence_limit), sigdigits=2)) V\n",
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
                "\t\tx: $(round(abs(refinement_limits[1]), sigdigits=2))\n",
                "\t\ty: $(round(abs(refinement_limits[2]), sigdigits=2))\n",
                "\t\tz: $(round(abs(refinement_limits[3]), sigdigits=2))\n",
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
                fssrb = PotentialSimulationSetupRB(detector, grid, potential, weighting_potential_contact_id = channel_id);
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

    return PotentialSimulationSetup{T, N, S}( Grid(fssrb), potential, pointtypes, ChargeDensityArray(fssrb), DielektrikumDistributionArray(fssrb)  )
end
