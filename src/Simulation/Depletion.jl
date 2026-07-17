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
        Umin::RealQuantity = minimum(broadcast(c -> c.potential, sim.detector.contacts)),
        Umax::RealQuantity = maximum(broadcast(c -> c.potential, sim.detector.contacts));
        contact_id::Int = determine_bias_voltage_contact_id(sim.detector),
        tolerance::RealQuantity = 0.1u"V",
        verbose::Bool = true,
        check_for_depletion::Bool = true) where {T <: AbstractFloat}

Estimates the potential needed to fully deplete the detector in a given [`Simulation`](@ref)
at the [`Contact`](@ref) with id `contact_id` by bisection method.
The default `contact_id` is determined automatically via `determine_bias_voltage_contact_id(sim.detector)`.
The default searching range of potentials is set by the extrema of contact potentials.

## Arguments 
* `sim::Simulation{T}`: [`Simulation`](@ref) for which the depletion voltage should be determined.
* `Umin::Real`: The minimum value of the searching range. If no units are given, this value is parsed in units of `$(internal_voltage_unit)`.
* `Umax::Real`: The maximum value of the searching range. If no units are given, this value is parsed in units of `$(internal_voltage_unit)`.
    
## Keywords
* `contact_id::Int`: The `id` of the [`Contact`](@ref) at which the potential is applied.
* `tolerance::Real`: The acceptable accuracy of results (default = 0.1V). If no units are given, this value is parsed in units of `$(internal_voltage_unit)`.
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
    U_min::RealQuantity = minimum(broadcast(c -> c.potential, sim.detector.contacts)),
    U_max::RealQuantity = maximum(broadcast(c -> c.potential, sim.detector.contacts));
    contact_id::Int = determine_bias_voltage_contact_id(sim.detector),
    tolerance::RealQuantity = 0.1u"V",
    verbose::Bool = true,
    check_for_depletion::Bool = true) where {T <: AbstractFloat}

    Umin::T = _parse_value(T, U_min, internal_voltage_unit) 
    Umax::T = _parse_value(T, U_max, internal_voltage_unit)
    tol::T  = _parse_value(T, tolerance, internal_voltage_unit)
    
    @assert !ismissing(sim.point_types) "Please calculate the electric potential first using `calculate_electric_potential!(sim)`"
    if check_for_depletion
        @assert is_depleted(sim.point_types) "This method only works for fully depleted simulations. Please increase the potentials in the configuration file to a greater value."
    end
    @assert Umax * Umin ≥ 0 "The voltage range needs to be positive or negative. Please adjust the voltage range."
    @assert abs(Umax - Umin) > tol "Umax - Umin < tolerance. Please change Umin, Umax or tolerance." 

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
    inside = findall((simDV.point_types .& 4 .> 0) .& (simDV.point_types .& inactive_layer_bit .== 0))
    U::T = NaN
    ϕ̃ = similar(ϕV)
    while Umax - Umin > tol
        U = (Umax + Umin)/2    
        ϕ̃ .= ϕρ .+ T(U) .* ϕV
        ϕmax = maximum(ϕ̃[inside])
        ϕmin = minimum(ϕ̃[inside])
        if U>0
            ϕmax - ϕmin < abs(U) ? Umax = U : Umin = U
        else
            ϕmax - ϕmin < abs(U) ? Umin = U : Umax = U
        end        
    end
    bulk = findall((sim.point_types.data .& bulk_bit .> 0) .& (sim.point_types.data .& inactive_layer_bit .== 0))
    U_min_max = filter(in(Urng), _find_depletion_voltage_candidates(ϕρ, ϕV, bulk))
    U = U > 0 ? max(U, U_min_max...) : min(U, U_min_max...)
    if verbose
        @info "The depletion voltage is around $(round(U, digits = Int(ceil(-log10(tol))))) ± $(tol) $(internal_voltage_unit) applied to contact $(contact_id)."
    end
    if (potential_range[2] - U) < tol || (U - potential_range[1]) < tol
        throw(ArgumentError("The depletion voltage ($(U * internal_voltage_unit)) is outside or too close to the edge of the search range $(potential_range .* internal_voltage_unit). Widen the range via `Umin`/`Umax`."))
    end
    return U * internal_voltage_unit
end


function _find_depletion_voltage_candidates(ϕᵨ::AbstractArray{T, 3}, ϕᵥ::AbstractArray{T, 3}, 
        bulk_points::Vector{CartesianIndex{3}}) where {T}

    L = length(bulk_points)
    Umin = zeros(L)
    Umax = zeros(L)
    indices = CartesianIndices(ϕᵥ)
    maxI = last(indices)
    minI = first(indices)
    Δ = oneunit(minI)
    for (i, x₀) in enumerate(bulk_points)
        ϕᵨ₀ = ϕᵨ[x₀]
        ϕᵥ₀ = ϕᵥ[x₀]
        neighbours = max(x₀ - Δ, minI):min(x₀ + Δ, maxI)
       _uminl = -Inf
       _uminr = Inf
       _umaxl = -Inf
       _umaxr = Inf
        for x in neighbours
            x₀ == x && continue
            δ1 = ϕᵨ₀ - ϕᵨ[x]
            δ2 = ϕᵥ₀ - ϕᵥ[x]
            isapprox(δ2, zero(T)) && continue
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
        Umin[i] = min((_uminl < _uminr) ? _uminl : Inf, (_umaxl < _umaxr) ? _umaxl : Inf)
        Umax[i] = max((_uminl < _uminr) ? _uminr : -Inf, (_umaxl < _umaxr) ? _umaxr : -Inf)
    end
    minimum(Umin), maximum(Umax)
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
