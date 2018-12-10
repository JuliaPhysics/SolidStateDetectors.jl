"""
    calculate_electric_potential(det::SolidStateDetector; <keyword arguments>)::Tuple{<:Grid, PointTypes}


Compute the electric potential for the given Detector `det` on an adaptive grid
through successive over relaxation.

There are serveral `<keyword arguments>` which can be used to tune the computation:
    
# Keywords
- `coordinates::Symbol`: the kind of the coordinate system of the grid. Right now only `:cylindrical` is possible.
- `convergence_limit::Real`: `convergence_limit` times the bias voltage sets the convergence limit of the relaxation. The convergence value is the absolute maximum difference of the potential between two iterations of all grid points. Default of `convergence_limit` is `5e-6` (times bias voltage).
- `max_refinements::Int`: number of maximum refinements. Default is `2`. Set it to `0` to switch off refinement.
- `refinement_limits::Vector{Real}`: vector of refinement limits for each dimension (in case of cylindrical coordinates the order is `r`, `φ`, `z`). A refinement limit (e.g. `refinement_limits[1]`) times the bias voltage of the detector `det` is the maximum allowed voltage difference between two neighbouring grid points in the respective dimension. When the difference is larger, new points are created inbetween. Default is `[1e-4, 1e-4, 1e-4]`.
- `min_grid_spacing::Vector{Real}`: vector of the mimimum allowed distance between two grid points for each dimension. For normal coordinates the unit is meter. For angular coordinates, the unit is radiance. It prevents the refinement to make the grid to fine. Default is [`1e-4`, `1e-2`, `1e-4`].
- `depletion_handling::Bool`: enables the handling of undepleted regions. Default is false.
- `nthreads::Int`: Number of threads to use in the computation. Default is `Base.Threads.nthreads()`. The environment variable `JULIA_NUM_THREADS` must be set appropriately before the Julia session was started (e.g. `export JULIA_NUM_THREADS=8` in case of bash).
- `sor_consts::Vector{<:Real}`: Two element array. First element contains the SOR constant for `r` = 0. Second contains the constant at the outer most grid point in `r`. A linear scaling is applied in between. First element should be smaller than the second one and both should be ∈ [1.0, 2.0]. Default is [1.4, 1.85].
- `return_2π_grid::Bool`: If set to `true`, the grid is extended to 2π in `φ` after it is computated. Default is `true`. This keyword is ignored if the simulation is in 2D. Use `extent_2D_grid_to_3D()` function to extend a 2D grid into 3D.
- `max_n_iterations::Int`: Set the maximum number of iterations which are performed after each grid refinement. Default is `10000`. If set to `-1` there will be no limit.
- `verbose::Bool=true`: Boolean whether info output is produced or not.

# Additional Information
- The function always turns the detector into a p-type detector to compute the potential. At the end it turns the signs to make it n-type again if it was n-type.
"""
function calculate_electric_potential(  det::SolidStateDetector;
                                        coordinates::Symbol=:cylindrical,
                                        convergence_limit::Real=5e-6,
                                        max_refinements::Int=3,
                                        refinement_limits::Vector{<:Real}=[1e-4, 1e-4, 1e-4],
                                        min_grid_spacing::Vector{<:Real}=[1e-4, 1e-2, 1e-4],
                                        init_grid_spacing::Real=2e-3, # 2 mm
                                        depletion_handling::Bool=false,
                                        nthreads::Int=Base.Threads.nthreads(),
                                        sor_consts::Vector{<:Real}=[1.4, 1.85],
                                        return_2π_grid::Bool=true,
                                        max_n_iterations::Int=10000,
                                        verbose::Bool=true,
                                        )::Tuple{<:Grid, PointTypes}
    T = get_precision_type(det)
    bias_voltage = T(det.segment_bias_voltages[1]);
    refinement_limits = T.(refinement_limits)
    min_grid_spacing = T.(min_grid_spacing)
    init_grid_spacing = T(init_grid_spacing)
    sor_consts = T.(sor_consts)
    refine::Bool = max_refinements > 0 ? true : false

    if nthreads > Base.Threads.nthreads()
        nthreads = Base.Threads.nthreads();
        @warn "`nthreads was set to `1`. The environment variable `JULIA_NUM_THREADS` must be set appropriately before the julia session is started."
    end

    cyclic = if 0 <= T(det.cyclic) <= T(2π)
        T(det.cyclic)
    else
        @warn "'cyclic=$(round(rad2deg(det.cyclic), digits = 0))°' is ∉ [0, 360] -> Set 'cyclic' to 360°."
        det.cyclic = T(2π)
        det.cyclic
    end
    only_2d::Bool = cyclic == 0 ? true : false

    n_φ_sym::Int = if !only_2d
        try
            Int(round(T(2π) / cyclic, digits = 3))
        catch err
            @warn "360° divided by 'cyclic=$(round(rad2deg(cyclic), digits = 0))°' does not give an integer -> Set 'cyclic' to 360°."
            1
        end
    else
        1
    end

    if n_φ_sym == 1 && !only_2d cyclic = T(2π) end
    n_φ_sym_info_txt = if n_φ_sym > 1
        "φ symmetry: cyclic = $(round(rad2deg(cyclic), digits = 0))° -> just calculating 1/$(n_φ_sym) of the detector in φ."
    else
        "φ symmetry: cyclic = $(round(rad2deg(cyclic), digits = 0))° -> no symmetry -> calculating for 360°."
    end
    if only_2d n_φ_sym_info_txt = "Detector is φ-symmetric -> 2D computation." end

    det_tmp = deepcopy(det)
    # always calculate for ptype detector -> if ntype change signs
    # segments are always set to 0 (normaly)
    if det_tmp.bulk_type == :ntype det_tmp.segment_bias_voltages[:] *= -1 end
    # For ptype detectors the (fixed space) charge carriere density is negative
    if det_tmp.charge_carrier_density_top > 0 det_tmp.charge_carrier_density_top *= -1 end
    if det_tmp.charge_carrier_density_bot > 0 det_tmp.charge_carrier_density_bot *= -1 end

    if bias_voltage == 0
        bias_voltage = if det.bulk_type == :ptype
            -abs(maximum( det_tmp.segment_bias_voltages[:]))
        else
            abs(minimum( det_tmp.segment_bias_voltages[:]))
        end
    end

    if verbose
        println("Electric Potential Calculation\n",
                "Bulk type: $(det.bulk_type)\n",
                "Bias voltage: $bias_voltage V\n",
                "$n_φ_sym_info_txt\n",
                "Precision: $T\n",
                "Convergence limit: $convergence_limit => $(round(abs(bias_voltage * convergence_limit), sigdigits=2)) V\n",
                "Threads: $nthreads\n",
                "Refine? -> $refine\n",
                "Refinement parameters:\n",
                "\tmaximum number of refinements: $(max_refinements)\n",
                "\tminimum grid spacing:\n",
                "\t\tr: $(min_grid_spacing[1]) m\n",
                "\t\tφ: $(min_grid_spacing[2]) rad\n",
                "\t\tz: $(min_grid_spacing[3]) m\n",
                "\tRefinement limits:\n",
                "\t\tr: $(refinement_limits[1]) -> $(round(abs(bias_voltage * refinement_limits[1]), sigdigits=2)) V\n",
                "\t\tφ: $(refinement_limits[2]) -> $(round(abs(bias_voltage * refinement_limits[2]), sigdigits=2)) V\n",
                "\t\tz: $(refinement_limits[3]) -> $(round(abs(bias_voltage * refinement_limits[3]), sigdigits=2)) V\n",
                ""
        )
    end


    cg = CylindricalGrid(det_tmp, init_grid_spacing=init_grid_spacing, cyclic=cyclic, only_2d = only_2d );
    pwrb = PrecalculatedWeightsCylindricalRedBlack(det_tmp, cg, sor_consts=sor_consts, only_2d = only_2d);
    rbgrid = CylindricalRedBlackGrid(cg, only_2d=only_2d);

    update_till_convergence!(rbgrid, pwrb, convergence_limit, bias_voltage, depletion_handling=depletion_handling, nthreads=nthreads, max_n_iterations=max_n_iterations)
    overwrite_potential!(cg, rbgrid);

    refinement_counter::Int = 0
    if refine
        new_inds = check_for_refinement(cg, pwrb, abs.(refinement_limits .* bias_voltage), min_grid_spacing, info_output=false)
        while !isempty(new_inds[1]) || !isempty(new_inds[2]) || !isempty(new_inds[3])
            refinement_counter += 1
            if refinement_counter <= max_refinements
                if verbose print("Refinement $refinement_counter:\t") end
                add_points_and_interpolate!(cg, new_inds...);
                if verbose println("New grid size: $(size(cg.potential))") end
                pwrb = PrecalculatedWeightsCylindricalRedBlack(det_tmp, cg, sor_consts=sor_consts, only_2d = only_2d);
                rbgrid = CylindricalRedBlackGrid(cg, only_2d=only_2d);
                update_till_convergence!(rbgrid, pwrb, convergence_limit, bias_voltage, depletion_handling=depletion_handling, nthreads=nthreads, max_n_iterations=max_n_iterations)
                overwrite_potential!(cg, rbgrid);
                new_inds = check_for_refinement(cg, pwrb, abs.(refinement_limits .* bias_voltage), min_grid_spacing, info_output=false)
            else
                break
            end
        end
    end

    pointypes = get_PointTypes(cg, pwrb)

    bias_voltage_abs::T = abs(bias_voltage)

    @inbounds for i in eachindex(pointypes)
        if (pointypes[i] & update_bit == 0)
            pointypes[i] = PointType(0)
        else
            if (pointypes[i] & undepleted_bit > 0) && (0 < abs(cg.potential[i]) < bias_voltage_abs)
                pointypes[i] += bubble_bit
            end
            if (pointypes[i] & pn_junction_bit == 0)
                if pointypes[i] & undepleted_bit > 0 pointypes[i] -= undepleted_bit end
                if pointypes[i] & bubble_bit > 0 pointypes[i] -= bubble_bit end
            end
        end
    end
    pts = PointTypes{T}(cg.r, cg.φ, cg.z, pointypes )

    # turn n-type detector to n-type again (since it was calculated as p-type)
    abs_max::T = abs(maximum(cg.potential))
    abs_min::T = abs(minimum(cg.potential))
    if det_tmp.bulk_type == :ntype cg.potential *= -1 end

    # Checks
    if bias_voltage_abs < abs_max + abs_min
    # if abs_max > abs(det.segment_bias_voltages[1])
        @warn "Detector not fully depleted at this bias voltage ($(det.segment_bias_voltages[1]) V). At least on grid point has a higher (or lower) potential value ($(sign(det.segment_bias_voltages[1]) * abs_max) V). This should not be."
    end
    if !in(det.segment_bias_voltages[1], cg.potential)
        @warn "Something is wrong. Applied bias voltage ($(det.segment_bias_voltages[1]) V) does not occur in the grid. This should not be since it is a fixed value."
    end


    rg = if !only_2d
        return_2π_grid ? extent_in_φ_n_times(cg, n_φ_sym) : cg
    else
        return cg, pts
    end
    return rg, pts
end
