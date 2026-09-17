"""
    set_point_type_depletion_handling!(sim::Simulation, depletion_handling::Bool)

Update the `depletion_handling` flag of the `PointTypes` saved in a simulation `sim` to the given value.

## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) for which the [`PointTypes`](@ref) depletion handling should be updated.
* `depletion_handling::Bool`: Value to which it should be updated to.
"""
function set_point_type_depletion_handling!(sim::Simulation, depletion_handling::Bool)
    sim.point_types = PointTypes(sim.point_types.data, sim.point_types.grid, depletion_handling)
    nothing
end


# """
#     _adapt_weighting_potential_to_electric_potential_grid!(sim::Simulation, contact_id::Int)
# 
# Interpolates the [`WeightingPotetial`](@ref) of the [`Contact`](@ref) with id `contact_id`
# onto the [`Grid`](@ref) of the [`ElectricPotential`](@ref) and updates it until it converges.
# 
# ## Arguments 
# * `sim::Simulation{T}`: [`Simulation`](@ref) for which the [`WeightingPotential`] should be adapted.
# * `contact_id::Int`: The `id` of the [`Contact`](@ref) at which the [`WeightingPotential`] should be adapted.
# """
function _adapt_weighting_potential_to_electric_potential_grid!(sim::Simulation, contact_id::Int)
    @assert !ismissing(sim.electric_potential) "Electric potential missing"
    if ismissing(sim.weighting_potentials[contact_id]) calculate_weighting_potential!(sim, contact_id, verbose = false) end
    if sim.weighting_potentials[contact_id].grid != sim.electric_potential.grid
        sim.weighting_potentials[contact_id] = sim.weighting_potentials[contact_id][sim.electric_potential.grid]
        #apply_initial_state!(sim, WeightingPotential, contact_id, sim.electric_potential.grid)
        update_till_convergence!(sim, WeightingPotential, contact_id)
    end
end

"""
    estimate_depletion_voltage(sim::Simulation{T},
        Umin::T = minimum(broadcast(c -> c.potential, sim.detector.contacts)),
        Umax::T = maximum(broadcast(c -> c.potential, sim.detector.contacts));
        contact_id::Int = determine_bias_voltage_contact_id(sim.detector),
        tolerance::AbstractFloat = 1e-1,
        verbose::Bool = true,
        check_for_depletion::Bool = true) where {T <: AbstractFloat}

Estimates the potential needed to fully deplete the detector in a given [`Simulation`](@ref)
at the [`Contact`](@ref) with id `contact_id` by bisection method.
The default `contact_id` is determined automatically via `determine_bias_voltage_contact_id(sim.detector)`.
The default searching range of potentials is set by the extrema of contact potentials.

## Arguments 
* `sim::Simulation{T}`: [`Simulation`](@ref) for which the depletion voltage should be determined.
* `Umin::T`: The minimum value of the searching range.
* `Umax::T`: The maximum value of the searching range.
    
## Keywords
* `contact_id::Int`: The `id` of the [`Contact`](@ref) at which the potential is applied.
* `tolerance::AbstractFloat`: The acceptable accuracy of results. Default is 1e-1.
* `verbose::Bool = true`: Activate or deactivate additional info output. Default is `true`.
* `check_for_depletion::Bool = true`: If `true` (default), assert via [`is_depleted`](@ref) that the
    simulation is fully depleted before estimating, since the method is only valid for a depleted detector.
    Set to `false` to skip this assertion.

## Example
```julia
using SolidStateDetectors
sim = Simulation(SSD_examples[:InvertedCoax])
calculate_electric_potential!(sim)
estimate_depletion_voltage(sim)
```

!!! note
    The accuracy of the result depends on the precision of the initial simulation.

!!! note
    This function performs two 2D or 3D field calculations, depending on `sim`.\\
    Thus, keep in mind that is might consume some memory. 
    
See also [`is_depleted`](@ref).

"""
function estimate_depletion_voltage(sim::Simulation{T},
    Umin::T = minimum(broadcast(c -> c.potential, sim.detector.contacts)),
    Umax::T = maximum(broadcast(c -> c.potential, sim.detector.contacts));
    contact_id::Int = determine_bias_voltage_contact_id(sim.detector),
    tolerance::AbstractFloat = 1e-1,
    verbose::Bool = true,
    check_for_depletion::Bool = true) where {T <: AbstractFloat}

    @assert !ismissing(sim.point_types) "Please calculate the electric potential first using `calculate_electric_potential!(sim)`"
    if check_for_depletion
        @assert is_depleted(sim.point_types) "This method only works for fully depleted simulations. Please increase the potentials in the configuration file to a greater value."
    end
    @assert Umax * Umin ≥ 0 "The voltage range needs to be positive or negative. Please adjust the voltage range."
    @assert abs(Umax - Umin) > tolerance "Umax - Umin < tolerance. Please change Umin, Umax or tolerance." 

    if all(sim.q_eff_imp.data .== 0) && all(sim.q_eff_fix.data .== 0)
        @warn "The detector seems to have no impurities. Therefore, the depletion voltage is 0 V"
        return zero(T) * u"V"
    end

    potential_range::Tuple{T, T} = (Umin, Umax)
    if verbose
        @info "Looking for the depletion voltage applied to contact $(contact_id) "*
        "in the range $(potential_range.*u"V")."
    end  

    simDV = deepcopy(sim)
    if ismissing(simDV.weighting_potentials[contact_id]) || simDV.weighting_potentials[contact_id].grid != simDV.electric_potential.grid
        _adapt_weighting_potential_to_electric_potential_grid!(simDV, contact_id)
    end

    # ϕρ is the electric potential resulting only from the impurity density
    ϕρ = simDV.electric_potential.data
    ϕV = simDV.weighting_potentials[contact_id].data   

    ϕρ .-= simDV.detector.contacts[contact_id].potential .* ϕV
    Urng = Umin..Umax

    bulk = findall((sim.point_types.data .& bulk_bit .> 0))
    U_cand = filter(in(Urng), _find_depletion_voltage_candidates(ϕρ, ϕV, bulk))
    U::T = if length(U_cand) == 0
        inside = findall(simDV.point_types .& 4 .> 0)
        U_bisection = _find_depletion_voltage_by_bisection(ϕρ, ϕV, inside, bulk, Umin, Umax, T(tolerance))
        if verbose
            @info """The local extremum scan returned no candidates (or no unique candidate). Using alternate bisection calculation.
            The depletion voltage is $(round(U_bisection, digits = Int(ceil(-log10(tolerance))))) ± $(tolerance) $(internal_voltage_unit) applied to contact $(contact_id)."""
        end
        U_bisection
    else
        U_cand = mean(Urng) > 0 ? max(U_cand...) : min(U_cand...)
        if verbose
            @info "The depletion voltage is $(round(U_cand, digits = 2)) $(internal_voltage_unit) applied to contact $(contact_id)."
        end
        U_cand
    end
    if (potential_range[2] - U) < tolerance || (U - potential_range[1]) < tolerance
        throw(ArgumentError("The depletion voltage ($(U * internal_voltage_unit)) is outside or too close to the edge of the search range $(potential_range .* internal_voltage_unit). Widen the range via `Umin`/`Umax`."))
    end
    return U * u"V"
end

function _has_local_extremum(ϕ::AbstractArray{T, 3}, bulk_points::Vector{CartesianIndex{3}}; verbose::Bool = false) where {T}
    indices = CartesianIndices(ϕ)
    maxI = last(indices)
    minI = first(indices)
    Δ = oneunit(minI)
    for x₀ in bulk_points
        is_local_minimum = true
        is_local_maximum = true
        ϕ₀ = ϕ[x₀]
        neighbours = max(x₀ - Δ, minI):min(x₀ + Δ, maxI)
        for x in neighbours
            Δx = sum(abs.((x₀ - x).I))
            ((x₀ == x) || Δx > 1) && continue
            is_local_minimum &= ϕ[x] > ϕ₀
            is_local_maximum &= ϕ[x] < ϕ₀
            (is_local_minimum || is_local_maximum) || break
        end
        if is_local_minimum || is_local_maximum
            verbose && @info "interior local $(is_local_minimum ? "minimum" : "maximum") of the potential at $(x₀.I), ϕ = $ϕ₀"
            return true
        end
    end
    return false
end

"""
    _find_depletion_voltage_candidates(ϕᵨ::AbstractArray{T, 3}, ϕᵥ::AbstractArray{T, 3},
        bulk_points::Vector{CartesianIndex{3}}) where {T}

Determines the range of bias voltages `U` for which at least one of the `bulk_points`
is a local extremum of the electric potential, i.e. for which the detector is not fully depleted.

The electric potential is linear in the bias voltage,

`ϕ(U) = ϕᵨ + U·ϕᵥ`,

so for a bulk point `x₀` and one of its grid neighbours `x` the sign of the difference
`ϕ(x₀) - ϕ(x) = δ1 + U·δ2`, with `δ1 = ϕᵨ[x₀] - ϕᵨ[x]` and `δ2 = ϕᵥ[x₀] - ϕᵥ[x]`,
flips exactly once at `U = -δ1/δ2`. Intersecting these thresholds over all neighbours in the
surrounding `3×3×3` block yields, for every bulk point, the interval of `U` in which it is a
local maximum (undepleted n-type region) and the interval in which it is a local minimum
(undepleted p-type region). Neighbours with `ϕᵥ[x] ≈ ϕᵥ[x₀]` have a `U`-independent ordering
and can rule out one of the two cases entirely.

The union of all valid intervals is the range of bias voltages for which the detector is undepleted.
Its extrema are returned: for a positive (negative) bias the detector is fully depleted for all
`U` above `Umax` (below `Umin`), so `Umax` (`Umin`) is the depletion voltage candidate.
[`estimate_depletion_voltage`](@ref) uses the candidate lying in its search range to refine the
result.

## Arguments
* `ϕᵨ::AbstractArray{T, 3}`: Electric potential resulting only from the impurity density, i.e. the
    full electric potential minus `V·ϕᵥ`, where `V` is the potential currently applied to the bias contact.
* `ϕᵥ::AbstractArray{T, 3}`: [`WeightingPotential`](@ref) of the bias contact, on the same grid as `ϕᵨ`.
* `bulk_points::Vector{CartesianIndex{3}}`: Indices of the grid points to check.

## Returns
A tuple `(Umin, Umax)` in units of `$(internal_voltage_unit)`. If no bulk point is a local extremum
for any `U`, the tuple `(Inf, -Inf)` is returned.
"""
function _find_depletion_voltage_candidates(ϕᵨ::AbstractArray{T, 3}, ϕᵥ::AbstractArray{T, 3}, 
        bulk_points::Vector{CartesianIndex{3}}) where {T}

    Umin_all = Inf
    Umax_all = -Inf
    indices = CartesianIndices(ϕᵥ)
    maxI = last(indices)
    minI = first(indices)
    Δ = oneunit(minI)
    for x₀ in bulk_points
        ϕᵨ₀ = ϕᵨ[x₀]
        ϕᵥ₀ = ϕᵥ[x₀]
        neighbours = max(x₀ - Δ, minI):min(x₀ + Δ, maxI)
       _uminl = -Inf # local min left bound, to be maximized for most stringent bound
       _uminr = Inf # local min right bound, to be minimized for most stringent bound
       _umaxl = -Inf # local max left bound, to be maximized for most stringent bound
       _umaxr = Inf # local max right bound, to be minimized for most stringent bound
        min_ok = true
        max_ok = true
        for x in neighbours
            x₀ == x && continue
            δ1 = ϕᵨ₀ - ϕᵨ[x]
            δ2 = ϕᵥ₀ - ϕᵥ[x]
            # min -> p-type: Find U range that solves δϕ = δ1 + U*δ2 < 0  
            # max -> n-type: Find U range that solves δϕ = δ1 + U*δ2 > 0 
            # TODO: Optimize the tolerance. Detectors with large regions of very small weighting potential, ϕᵥ, can return no candidates.
            if isapprox(ϕᵥ₀, ϕᵥ[x], atol = 20*eps(T)) # ordering independent of U
                δ1 >= 0 && (min_ok = false) # no local minimum for any U, local maximum for every U -> keep [_umaxl, _uminr] = [-Inf, Inf]
                δ1 <= 0 && (max_ok = false) # no local maximum for any U, local minimum for every U -> keep [_uminl, _uminr] = [-Inf, Inf]
                continue
            end
            U = -(δ1 / δ2) 
            if δ2 > 0
                _umaxl = max(_umaxl, U)
                _uminr = min(_uminr, U)
            end
            if δ2 < 0
                _umaxr = min(_umaxr, U)
                _uminl = max(_uminl, U)
            end
        end 
        max_ok &= _umaxl < _umaxr # if out of order -> no local maximum for any U
        min_ok &= _uminl < _uminr # if out of order -> no local minimum for any U
        # find greatest (lowest) U for which any bulk point is a local minimum (maximum) -> Above (below) this U depletion is achived
        for (ok, l, r) in ((max_ok, _umaxl, _umaxr), (min_ok, _uminl, _uminr))
            ok || continue  # only check valid bounds
            Umax_all = max(Umax_all, r)
            Umin_all = min(Umin_all, l)
        end
    end
    Umin_all, Umax_all
end

function _find_depletion_voltage_by_bisection(ϕᵨ::AbstractArray{T, 3}, ϕᵥ::AbstractArray{T, 3}, 
        inside_points::Vector{CartesianIndex{3}}, bulk_points::Vector{CartesianIndex{3}}, Umin::T, Umax::T, tol::T) where {T}
    U::T = NaN
    ϕ̃ = similar(ϕᵥ)
    while Umax - Umin > tol
        U = (Umax + Umin)/2    
        ϕ̃ .= ϕᵨ .+ T(U) .* ϕᵥ
        ϕmax = maximum(ϕ̃[inside_points])
        ϕmin = minimum(ϕ̃[inside_points])
        if U>0
            ϕmax - ϕmin < abs(U) && !_has_local_extremum(ϕ̃, bulk_points) ? Umax = U : Umin = U
        else
            ϕmax - ϕmin < abs(U) && !_has_local_extremum(ϕ̃, bulk_points) ? Umin = U : Umax = U
        end        
    end
    U
end

"""
    adjust_impurity_and_electric_potential_to_match_depletion!(sim::Simulation{T}, dep::RealQuantity,
        U_min::RealQuantity = minimum(broadcast(c -> c.potential, sim.detector.contacts)),
        U_max::RealQuantity = maximum(broadcast(c -> c.potential, sim.detector.contacts));
        contact_id::Int = determine_bias_voltage_contact_id(sim.detector),
        verbose::Bool = true,
        reconverge_electric_potential::Bool = true,
        kwargs...) where {T <: AbstractFloat}

Rescales the impurity density of a [`Simulation`](@ref) in place so that its depletion voltage
matches the target `dep`, **without a full electric potential re-solve** and **without changing the bias**.

The `impurity_density_model` is scaled by `f = dep / dep_sim`, where `dep_sim` is the current
depletion voltage estimated via [`estimate_depletion_voltage`](@ref). Since the electric potential
is linear in the impurity density, the stored `sim.electric_potential` is updated by scaling its
impurity contribution by `f` (the bias contribution, given by the [`WeightingPotential`](@ref) of
the bias contact, is left unchanged). With the bias voltage `V` held fixed, the stored potential is
updated as

`ϕ → f·ϕ + (1 - f)·V·ϕV`,    where    `f = dep / dep_sim`

and `ϕV` is the weighting potential of the bias contact.

!!! note
    The [`WeightingPotential`](@ref) of the bias contact (`contact_id`) may be modified by this
    function: if it is missing or does not share the grid of the [`ElectricPotential`](@ref), it
    is calculated/mapped onto that grid.

!!! note
    This modifies `sim.electric_potential`. If `sim.electric_field` has already been calculated, it
    must be recomputed via [`calculate_electric_field!`](@ref) to reflect the change.

The target `dep` must share the same (non-zero) sign as the current bias voltage and be smaller
in magnitude (i.e. the current bias over-depletes the detector at the target depletion voltage),
otherwise an `ArgumentError` is thrown.

## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) to be adapted in place.
* `dep::RealQuantity`: Target depletion voltage. If no units are given, this value is parsed in units of `$(internal_voltage_unit)`.
* `U_min::RealQuantity`, `U_max::RealQuantity`: Search range forwarded to [`estimate_depletion_voltage`](@ref)
    when determining the current depletion voltage `dep_sim`. Defaults to the extrema of the contact potentials.

## Keywords
* `contact_id::Int`: The `id` of the [`Contact`](@ref) at which the bias voltage is applied.
    The default is determined automatically via `determine_bias_voltage_contact_id(sim.detector)`.
* `verbose::Bool = true`: Activate or deactivate additional info output. Default is `true`.
* `reconverge_electric_potential::Bool = true`: If `true` (default), the [`ElectricPotential`](@ref) is
    relaxed to convergence via `update_till_convergence!` after the analytical superposition (using the
    superposed field as the initial state).
    Set to `false` to keep the purely analytical superposition result.

Additional `kwargs...` are passed on to [`estimate_depletion_voltage`](@ref).

!!! note
    When `reconverge_electric_potential = true` (the default), the field is re-relaxed by the SOR
    solver, so a subsequent call to [`estimate_depletion_voltage`](@ref) may return a value slightly
    different (usually a few volts) from the target `dep`. Set it to `false` for the exact analytical match.

Returns the impurity scaling factor `f = dep / dep_sim` that was applied.

See also [`adjust_bias_and_electric_potential!`](@ref).
"""
function adjust_impurity_and_electric_potential_to_match_depletion!(sim::Simulation{T}, dep::RealQuantity,
    U_min::RealQuantity = minimum(broadcast(c -> c.potential, sim.detector.contacts)),
    U_max::RealQuantity = maximum(broadcast(c -> c.potential, sim.detector.contacts));
    contact_id::Int = determine_bias_voltage_contact_id(sim.detector),
    verbose::Bool = true,
    reconverge_electric_potential::Bool = true,
    kwargs...) where {T <: AbstractFloat}

    dep::T = _parse_value(T, dep, internal_voltage_unit)
    V = sim.detector.contacts[contact_id].potential
    dep * V <= 0 && throw(ArgumentError("The depletion voltage ($(dep)$(internal_voltage_unit)) and operating voltage ($(V)$(internal_voltage_unit)) must have the same (non-zero) sign."))
    abs(V) <= abs(dep) && throw(ArgumentError("The operating voltage ($(V)$(internal_voltage_unit)) must exceed the depletion voltage ($(dep)$(internal_voltage_unit)) in magnitude. The detector must be over-depleted to superpose the contributions of the bias voltage and impurity density to the electric potential."))

    if ismissing(sim.weighting_potentials[contact_id]) || sim.weighting_potentials[contact_id].grid != sim.electric_potential.grid
        _adapt_weighting_potential_to_electric_potential_grid!(sim, contact_id)
    end
    dep_sim = _parse_value(T, estimate_depletion_voltage(sim, U_min, U_max; contact_id, verbose, kwargs...), internal_voltage_unit)
    f = T(dep / dep_sim)
    if verbose
        @info "Scaling `impurity_density_model` by $f so that the depletion voltage matches $(dep)$(internal_voltage_unit)."
    end
    ϕV = sim.weighting_potentials[contact_id].data
    # ϕ = f·ϕρ + V·ϕV = f·(ϕ - V·ϕV) + V·ϕV = f·ϕ + V·(1 - f)·ϕV
    sim.electric_potential.data .= f .* sim.electric_potential.data .+ (V * (1 - f)) .* ϕV
    sim.detector = SolidStateDetector(sim.detector, f * sim.detector.semiconductor.impurity_density_model)
    if reconverge_electric_potential
        update_till_convergence!(sim, ElectricPotential; verbose)
        mark_bits!(sim)
    end
    if verbose && !ismissing(sim.electric_field)
        @info "The impurity density was scaled to match the depletion and the electric potential was modified accordingly to reflect this change.\n" *
              "Please run `calculate_electric_field!` to reflect this change also in the electric field."
    end
    f
end

"""
    adjust_bias_and_electric_potential!(sim::Simulation{T}, bias::RealQuantity;
        contact_id::Int = determine_bias_voltage_contact_id(sim.detector),
        verbose::Bool = true,
        check_against_depletion_voltage::Bool = true,
        reconverge_electric_potential::Bool = true,
        kwargs...) where {T <: AbstractFloat}

Swaps the bias contact potential of a [`Simulation`](@ref) to `bias` in place, 
**without a full electric potential re-solve** and **without changing the impurity density**.

Since the electric potential is linear in the applied bias, the stored potential is updated as

`ϕ → ϕ + (bias - V_old)·ϕV`

where `V_old` is the current contact potential and `ϕV` is the [`WeightingPotential`](@ref) of the
bias contact.

!!! note
    The [`WeightingPotential`](@ref) of the bias contact (`contact_id`) may be modified by this
    function: if it is missing or does not share the grid of the [`ElectricPotential`](@ref), it
    is calculated/mapped onto that grid.

!!! note
    This modifies `sim.electric_potential`. If `sim.electric_field` has already been calculated, it
    must be recomputed via [`calculate_electric_field!`](@ref) to reflect the change.

## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) to be adapted in place.
* `bias::RealQuantity`: Target bias voltage. If no units are given, this value is parsed in units of `$(internal_voltage_unit)`.

## Keywords
* `contact_id::Int`: The `id` of the [`Contact`](@ref) at which the bias voltage is applied.
    The default is determined automatically via `determine_bias_voltage_contact_id(sim.detector)`.
* `verbose::Bool = true`: Activate or deactivate additional info output. Default is `true`.
* `check_against_depletion_voltage::Bool = true`: If `true`, the depletion voltage is estimated via
    [`estimate_depletion_voltage`](@ref) and an `ArgumentError` is thrown unless `bias` shares its
    (non-zero) sign and exceeds it in magnitude (i.e. the detector is fully depleted at `bias`). Set
    to `false` to skip this check, e.g. when the depletion voltage has just been set via
    [`adjust_impurity_and_electric_potential_to_match_depletion!`](@ref).
* `reconverge_electric_potential::Bool = true`: If `true` (default), the [`ElectricPotential`](@ref) is
    relaxed to convergence via `update_till_convergence!` after the analytical superposition (using the
    superposed field as the initial state).
    Set to `false` to keep the purely analytical superposition result.

Additional `kwargs...` are passed on to [`estimate_depletion_voltage`](@ref).

See also [`adjust_impurity_and_electric_potential_to_match_depletion!`](@ref).
"""
function adjust_bias_and_electric_potential!(sim::Simulation{T}, bias::RealQuantity;
    contact_id::Int = determine_bias_voltage_contact_id(sim.detector),
    verbose::Bool = true,
    check_against_depletion_voltage::Bool = true,
    reconverge_electric_potential::Bool = true,
    kwargs...) where {T <: AbstractFloat}

    V = sim.detector.contacts[contact_id].potential
    bias = _parse_value(T, bias, internal_voltage_unit)
    if check_against_depletion_voltage 
        dep = _parse_value(T, estimate_depletion_voltage(sim; contact_id, verbose, kwargs...), internal_voltage_unit) 
        dep * bias <= 0 && throw(ArgumentError("The depletion voltage ($(dep)$(internal_voltage_unit)) and operating voltage ($(bias)$(internal_voltage_unit)) must have the same (non-zero) sign."))
        abs(bias) <= abs(dep) && throw(ArgumentError("The operating voltage ($(bias)$(internal_voltage_unit)) must exceed the depletion voltage ($(dep)$(internal_voltage_unit)) in magnitude. The detector must be over-depleted to superpose the contributions of the bias voltage and impurity density to the electric potential."))
    end

    if ismissing(sim.weighting_potentials[contact_id]) || sim.weighting_potentials[contact_id].grid != sim.electric_potential.grid
        _adapt_weighting_potential_to_electric_potential_grid!(sim, contact_id)
    end
    V = sim.detector.contacts[contact_id].potential
    if verbose
        @info "Shifting the bias on `sim.detector.contacts[$contact_id]` from $(V)$(internal_voltage_unit) to $(bias)$(internal_voltage_unit)."
    end
    ϕV = sim.weighting_potentials[contact_id].data
    sim.electric_potential.data .+= (bias - V) .* ϕV
    sim.detector = SolidStateDetector(sim.detector, contact_id = contact_id, contact_potential = bias)
    if reconverge_electric_potential
        update_till_convergence!(sim, ElectricPotential; verbose)
        mark_bits!(sim)
    end
    if verbose && !ismissing(sim.electric_field)
        @info "The bias voltage was adjusted and the electric potential was modified accordingly to reflect this change.\n" *
              "Please run `calculate_electric_field!` to reflect this change also in the electric field."
    end
    nothing
end

#=
"""
    old_estimate_depletion_voltage( sim::Simulation{T}, contact_id::Int, field_sim_settings = (verbose = true,))::T

Estimates the full depletion voltage, U_D, of a detector in a given [`Simulation`](@ref).\\
That is is the voltage to fully deplete the detector.\\
This is done by calculating two electric potentials with different boundary conditions:\\
    1) Only the electric potential coming from the impurity density alone: `EP_i`\\
    2) The electric potential without an impurity density: `EP_0`\\
Then, the superpostion `EP_i + U_D * EP_0` is iteratively solved over all grid points\\
to determine `U_D` which is the voltage where the gradient of the suposition field, the electric field,\\
becomes nowhere 0 anywhere inside the semiconductor. 

## Arguments 
* `sim::Simulation{T}`: [`Simulation`](@ref) of the detector for which the depletion voltage should be determined.

## Keywords
* `bias_voltage_contact_id::Int`: The `id` of the [`Contact`](@ref) at which the potential is applied.\\
As default it is tried to determine this automatically via `determine_bias_voltage_contact_id(sim.detector)`.
* `field_sim_settings::NamedTuple = (verbose = false,)`: NamedTuple of simulation settings merged with the\\
default settings (listed underneath) and passed further to the field calculation functions.

The following default settings are used for the field simulations performed in this function:
```julia
(
    convergence_limit = 1e-7,
    max_tick_distance = 2.0u"mm",
    refinement_limits = [0.2, 0.1, 0.05, 0.025, 0.01], 
    sor_consts = (1.0, 1.0),
    use_nthreads = max_threads > 16 ? 16 : max_threads,
    n_iterations_between_checks = 20, 
    depletion_handling = false
)
```

## Example 
```julia  
using SolidStateDetectors
sim = Simulation(SSD_examples[:InvertedCoax])
old_estimate_depletion_voltage(sim, field_sim_settings = (verbose = true,))
```

!!! note
    The accuracy of the result depends on the precision of the initial simulation.

!!! note
    This function performs two 2D or 3D, depending on `sim`, field calculations.\\
    Thus, keep in mind that is might consume some memory. 
    
See also [`is_depleted`](@ref).
"""

function old_estimate_depletion_voltage(
        sim::Simulation{T,CS}, 
        field_sim_settings::NamedTuple = (verbose = false,);
        bias_voltage_contact_id::Int = determine_bias_voltage_contact_id(sim.detector)
    ) where {T <: AbstractFloat, CS}
    
    simDV = Simulation{T,CS}(
        deepcopy(sim.config_dict),
        deepcopy(sim.input_units),
        deepcopy(sim.medium),
        deepcopy(sim.detector),
        deepcopy(sim.world),
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        [missing for c in sim.weighting_potentials],
        missing
    )

    simDV.detector = SolidStateDetector(simDV.detector, contact_id = bias_voltage_contact_id, contact_potential = 0)

    max_threads = Base.Threads.nthreads()

    fss = merge((
        convergence_limit = 1e-7,
        max_tick_distance = max_tick_distance_default(simDV.world),
        refinement_limits = [0.2, 0.1, 0.05, 0.025, 0.01], 
        sor_consts = (1.2, 1.6),
        use_nthreads = max_threads > 16 ? 16 : max_threads,
        n_iterations_between_checks = 100, 
        depletion_handling = false
    ), field_sim_settings)

    calculate_electric_potential!(simDV; fss...)
    calculate_weighting_potential!(simDV, bias_voltage_contact_id; fss...)
    _adapt_weighting_potential_to_electric_potential_grid!(simDV, bias_voltage_contact_id)
    
    old_estimate_depletion_voltage(
        simDV.weighting_potentials[bias_voltage_contact_id],
        simDV.electric_potential,
        simDV.point_types, 
    )
end
        
function old_estimate_depletion_voltage(
    zero_imp_ep::ScalarPotential{T}, 
    only_imp_ep::ScalarPotential{T}, 
    point_types::PointTypes
    ) where {T}
    depletion_voltage_field = zeros(T, size(point_types.grid))
    
    _estimate_depletion_voltage_factor!(
        depletion_voltage_field, 
        point_types,
        only_imp_ep.data,
        zero_imp_ep.data
    )
    min, max = extrema(depletion_voltage_field)
    abs_depl_volt = (max - min) * u"V" # full depletion voltage
    sign_depl_volt = abs(min) > abs(max) ? -1 : 1
    sign_depl_volt * abs_depl_volt
end


function _estimate_depletion_voltage_factor(
        center_only_imp_value::T, 
        center_zero_imp_value::T, 
        neighbor_only_imp_values::AbstractVector{T}, 
        neighbor_zero_imp_values::AbstractVector{T}
    ) where {T}
    estimation = -center_only_imp_value / center_zero_imp_value 
    neighbor_values = neighbor_zero_imp_values .* estimation .+ neighbor_only_imp_values
    vmin, vmax = extrema(neighbor_values)
    center_value = center_zero_imp_value * estimation + center_only_imp_value
    if center_value < vmin
        estimation += (vmin - center_value) / center_zero_imp_value 
    elseif center_value > vmax
        estimation -= (vmax - center_value) / center_zero_imp_value
    end
    # neighbor_values = neighbor_zero_imp_values .* estimation .+ neighbor_only_imp_values
    # center_value = center_zero_imp_value * estimation + center_only_imp_value
    # vmin, vmax = extrema(neighbor_values)
    # @assert (center_value < vmin || center_value > vmax) 
    return estimation
end

function _replace_NaN_with_minimum(v::AbstractVector)
    min = minimum(filter(x -> !isnan(x), v))
    map(x -> isnan(x) ? min : x, v)
end

function _estimate_depletion_voltage_factor!(
        depletion_voltage_field::AbstractArray{T}, 
        point_types, 
        only_imp_ep, 
        zero_imp_ep
    ) where {T}
    @inbounds begin
        grid_size = size(point_types)
        for i3 in 1:grid_size[3]
            for i2 in 1:grid_size[2]
                for i1 in 1:grid_size[1]
                    if is_pn_junction_point_type(point_types[i1, i2, i3])
                        center_only_imp_value = only_imp_ep[i1, i2, i3]
                        center_zero_imp_value = zero_imp_ep[i1, i2, i3]
                        neighbor_only_imp_values::SVector{6, T} = _replace_NaN_with_minimum(@SVector T[
                            i1 < grid_size[1] && is_pn_junction_point_type(point_types[i1+1, i2, i3]) ? only_imp_ep[i1+1, i2, i3] : T(NaN),
                            i1 > 1 && is_pn_junction_point_type(point_types[i1-1, i2, i3]) ? only_imp_ep[i1-1, i2, i3] : T(NaN),
                            i2 < grid_size[2] && is_pn_junction_point_type(point_types[i1, i2+1, i3]) ? only_imp_ep[i1, i2+1, i3] : T(NaN),
                            i2 > 1 && is_pn_junction_point_type(point_types[i1, i2-1, i3]) ? only_imp_ep[i1, i2-1, i3] : T(NaN),
                            i3 < grid_size[3] && is_pn_junction_point_type(point_types[i1, i2, i3+1]) ? only_imp_ep[i1, i2, i3+1] : T(NaN),
                            i3 > 1 && is_pn_junction_point_type(point_types[i1, i2, i3-1]) ? only_imp_ep[i1, i2, i3-1] : T(NaN)
                        ])
                        neighbor_zero_imp_values::SVector{6, T} = _replace_NaN_with_minimum(@SVector T[
                            i1 < grid_size[1] && is_pn_junction_point_type(point_types[i1+1, i2, i3]) ? zero_imp_ep[i1+1, i2, i3] : T(NaN),
                            i1 > 1 && is_pn_junction_point_type(point_types[i1-1, i2, i3]) ? zero_imp_ep[i1-1, i2, i3] : T(NaN),
                            i2 < grid_size[2] && is_pn_junction_point_type(point_types[i1, i2+1, i3]) ? zero_imp_ep[i1, i2+1, i3] : T(NaN),
                            i2 > 1 && is_pn_junction_point_type(point_types[i1, i2-1, i3]) ? zero_imp_ep[i1, i2-1, i3] : T(NaN),
                            i3 < grid_size[3] && is_pn_junction_point_type(point_types[i1, i2, i3+1]) ? zero_imp_ep[i1, i2, i3+1] : T(NaN),
                            i3 > 1 && is_pn_junction_point_type(point_types[i1, i2, i3-1]) ? zero_imp_ep[i1, i2, i3-1] : T(NaN)
                        ])

                        depletion_voltage_field[i1, i2, i3] = _estimate_depletion_voltage_factor(
                            center_only_imp_value, center_zero_imp_value, neighbor_only_imp_values, neighbor_zero_imp_values
                        )
                    end
                end
            end
        end
    end
end
=#