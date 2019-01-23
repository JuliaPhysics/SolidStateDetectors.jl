"""
    calculate_weighting_potential(det::SolidStateDetector{T}; <keyword arguments>) where {T}


Compute the weighting potential for the given Detector `det` on an adaptive grid
through successive over relaxation.

There are serveral `<keyword arguments>` which can be used to tune the computation:
    
# Keywords
- `coordinates::Symbol`: the kind of the coordinate system of the grid. Right now only `:cylindrical` is possible.
- `convergence_limit::Real`: `convergence_limit` times the bias voltage sets the convergence limit of the relaxation. The convergence value is the absolute maximum difference of the potential between two iterations of all grid points. Default of `convergence_limit` is `5e-6` (times bias voltage).
- `max_refinements::Int`: number of maximum refinements. Default is `2`. Set it to `0` to switch off refinement.
- `refinement_limits::Vector{Real}`: vector of refinement limits for each dimension (in case of cylindrical coordinates the order is `r`, `θ`, `z`). A refinement limit (e.g. `refinement_limits[1]`) times the bias voltage of the detector `det` is the maximum allowed voltage difference between two neighbouring grid points in the respective dimension. When the difference is larger, new points are created inbetween. Default is `[1e-4, 1e-4, 1e-4]`.
- `min_grid_spacing::Vector{Real}`: vector of the mimimum allowed distance between two grid points for each dimension. For normal coordinates the unit is meter. For angular coordinates, the unit is radiance. It prevents the refinement to make the grid to fine. Default is [`1e-4`, `1e-2`, `1e-4`].
- `depletion_handling::Bool`: enables the handling of undepleted regions. Default is false.
- `use_nthreads::Int`: Number of threads to use in the computation. Default is `Base.Threads.nthreads()`. The environment variable `JULIA_NUM_THREADS` must be set appropriately before the Julia session was started (e.g. `export JULIA_NUM_THREADS=8` in case of bash).
- `sor_consts::Vector{<:Real}`: Two element array. First element contains the SOR constant for `r` = 0. Second contains the constant at the outer most grid point in `r`. A linear scaling is applied in between. First element should be smaller than the second one and both should be ∈ [1.0, 2.0]. Default is [1.4, 1.85].
- `return_2π_grid::Bool`: If set to `true`, the grid is extended to 2π in `θ` after it is computated. Default is `true`. This keyword is ignored if the simulation is in 2D. Use `extent_2D_grid_to_3D()` function to extend a 2D grid into 3D.
- `max_n_iterations::Int`: Set the maximum number of iterations which are performed after each grid refinement. Default is `10000`. If set to `-1` there will be no limit.
- `verbose::Bool=true`: Boolean whether info output is produced or not.
- `init_grid_spacing::Vector{<:Real}`: Initial spacing of the grid. Default is [2e-3, 2π / 72, 2e-3] <=> [2mm, 5 degree, 2mm ]

# Additional Information
- The function always turns the detector into a p-type detector to compute the potential. At the end it turns the signs to make it n-type again if it was n-type.
"""
function calculate_weighting_potential( detector::SolidStateDetector{T}, channel::Int;
                                        init_grid_spacing::Vector{<:Real} = [0.005, 20.0, 0.005], # 2mm, 10 degree, 2 mm
                                        grid::Grid{T, N, S} = Grid(detector, init_grid_spacing = init_grid_spacing, for_weighting_potential = true),
                                        convergence_limit::Real = 5e-6,
                                        max_refinements::Int = 3,
                                        refinement_limits::Vector{<:Real} = [1e-4, 1e-4, 1e-4], 
                                        min_grid_spacing::Vector{<:Real} = [1e-4, 1e-2, 1e-4],  # mm, degree, mm
                                        depletion_handling::Bool = false,
                                        use_nthreads::Int = Base.Threads.nthreads(),
                                        sor_consts::Vector{<:Real}=[1.4, 1.85],
                                        return_2π_grid::Bool=true,
                                        extend_2d_grid_to_n_points_in_θ::Int = 36,
                                        max_n_iterations::Int=10000,
                                        verbose::Bool=true,
                                        ) where {T, N, S}
    refinement_limits::Vector{T} = T.(refinement_limits)
    min_grid_spacing::Vector{T} = T.(min_grid_spacing)
    sor_consts::Vector{T} = T.(sor_consts)
    refine::Bool = max_refinements > 0 ? true : false
    only_2d::Bool = length(grid.axes[2]) == 1 ? true : false
    check_grid(grid)
    cyclic::T = grid.axes[2].interval.right
    n_θ_sym::Int = only_2d ? extend_2d_grid_to_n_points_in_θ : Int(round(T(2π) / cyclic, digits = 3)) 
    if use_nthreads > Base.Threads.nthreads()
        use_nthreads = Base.Threads.nthreads();
        @warn "`use_nthreads` was set to `1`. The environment variable `JULIA_NUM_THREADS` must be set appropriately before the julia session is started."
    end
    n_θ_sym_info_txt = if only_2d  
        "θ symmetry: Detector is θ-symmetric -> 2D computation."
    elseif n_θ_sym > 1 
        "θ symmetry: cyclic = $(round(rad2deg(cyclic), digits = 0))° -> calculating just 1/$(n_θ_sym) in θ of the detector."
    else
        "θ symmetry: cyclic = $(round(rad2deg(cyclic), digits = 0))° -> no symmetry -> calculating for 360°."
    end
    
    fssrb::PotentialSimulationSetupRB{T, 3, 4, :Cylindrical} = PotentialSimulationSetupRB(detector, grid, weighting_potential_channel_idx = channel);
    
    if verbose
        println("Weighting Potential Calculation\n",
                "Bulk type: $(detector.bulk_type)\n",
                "$n_θ_sym_info_txt\n",
                "Precision: $T\n",
                "Convergence limit: $convergence_limit => $(round(abs(convergence_limit), sigdigits=2)) V\n",
                "Threads: $use_nthreads\n",
                "Initial grid dimension: $(size(grid))\n",
                "Refine? -> $refine\n",
                "Refinement parameters:\n",
                "\tmaximum number of refinements: $(max_refinements)\n",
                "\tminimum grid spacing:\n",
                "\t\tr: $(min_grid_spacing[1]) m\n",
                "\t\tθ: $(min_grid_spacing[2]) rad\n",
                "\t\tz: $(min_grid_spacing[3]) m\n",
                "\tRefinement limits:\n",
                "\t\tr: $(refinement_limits[1]) -> $(round(abs(refinement_limits[1]), sigdigits=2)) V\n",
                "\t\tθ: $(refinement_limits[2]) -> $(round(abs(refinement_limits[2]), sigdigits=2)) V\n",
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
                fssrb = PotentialSimulationSetupRB(detector, grid, potential, weighting_potential_channel_idx = channel)
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
    
    pointtypes::Array{PointType, 3} = PointTypeArray(fssrb)
    # min_or_max_potential_for_bubble_marking::T = fssrb.bulk_is_ptype ? fssrb.minimum_applied_potential + 10 : fssrb.maximum_applied_potential - 10
    @inbounds for i in eachindex(pointtypes)
        if (pointtypes[i] & update_bit == 0)
            pointtypes[i] = PointType(0)
        else
            if (pointtypes[i] & undepleted_bit > 0) 
                # if fssrb.bulk_is_ptype 
                #     if (min_or_max_potential_for_bubble_marking < potential[i] < fssrb.maximum_applied_potential)
                #         pointtypes[i] += bubble_bit
                #     end
                # else
                #     if (fssrb.minimum_applied_potential < potential[i] < min_or_max_potential_for_bubble_marking)
                #         pointtypes[i] += bubble_bit
                #     end
                # end
            end
            if (pointtypes[i] & pn_junction_bit == 0)
                if pointtypes[i] & undepleted_bit > 0 pointtypes[i] -= undepleted_bit end
                # if pointtypes[i] & bubble_bit > 0 pointtypes[i] -= bubble_bit end
            end
        end
    end
        
    # rg = if !only_2d
        #     return_2π_grid ? extent_in_θ_n_times(cg, n_θ_sym) : cg
        # else
        #     return cg, pts
        # end
        # return rg, pts
        
    return PotentialSimulationSetup{T, N, S}( Grid(fssrb), potential, pointtypes, ChargeDensityArray(fssrb), DielektrikumDistributionArray(fssrb)  )
end
