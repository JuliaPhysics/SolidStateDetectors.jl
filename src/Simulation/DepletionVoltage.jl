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
    verbose::Bool = true)::T where {T <: AbstractFloat}

Estimates the potential (in V) needed to fully deplete the detector in a given [`Simulation`](@ref)
at the [`Contact`](@ref) with id `contact_id` by bisection method.
The default `contact_id` is determined automatically via `determine_bias_voltage_contact_id(sim.detector)`.
The default searching range of potentials is set by the extrema of contact potentials.

## Arguments 
* `sim::Simulation{T}`: [`Simulation`](@ref) for which the depletion voltage should be determined.
* `Umin::Int`: The minimum value of the searching range.
* `Umax::Int`: The maximum value of the searching range.
    
## Keywords
* `contact_id::Int`: The `id` of the [`Contact`](@ref) at which the potential is applied.
* `tolerance::AbstractFloat`: The acceptable accuracy of results. Default is 1e-1.
* `verbose::Bool = true`: Activate or deactivate additional info output. Default is `true`.

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
    This function performs two 2D or 3D, depending on `sim`, field calculations.\\
    Thus, keep in mind that is might consume some memory. 
    
See also [`is_depleted`](@ref).

"""
function estimate_depletion_voltage(sim::Simulation{T},
    Umin::T = minimum(broadcast(c -> c.potential, sim.detector.contacts)),
    Umax::T = maximum(broadcast(c -> c.potential, sim.detector.contacts));
    contact_id::Int = determine_bias_voltage_contact_id(sim.detector),
    tolerance::AbstractFloat = 1e-1,
    verbose::Bool = true)::T where {T <: AbstractFloat}

    potential_range::Tuple{T, T} = (Umin,Umax)

    @assert is_depleted(sim.point_types) "This method only works for fully depleted simulations. Please increase the potentials in the configuration file to a greater value."

    @assert Umax-Umin>tolerance "Umax - Umin < tolerance. Please change Umin, Umax or tolerance."

    simDV = deepcopy(sim)

    if verbose
        @info "Looking for the depletion voltage applied to contact $(contact_id) "*
        "in the range $(potential_range.*u"V")."
    end  

    if ismissing(simDV.weighting_potentials[contact_id]) || simDV.weighting_potentials[contact_id].grid != simDV.electric_potential.grid
        _adapt_weighting_potential_to_electric_potential_grid!(simDV, contact_id)
    end

    # ϕρ is the electric potential resulting only from the impurity density
    ϕρ = simDV.electric_potential.data
    ϕV = simDV.weighting_potentials[contact_id].data   

    ϕρ .-= simDV.detector.contacts[contact_id].potential .* ϕV

    inside = findall(simDV.point_types .& 4 .> 0)
    U::T = NaN

    while Umax-Umin>tolerance
        U = (Umax + Umin)/2    
        if U>0
            maximum(ϕρ[inside] .+ T(U) .* ϕV[inside]) - minimum(ϕρ[inside] .+ T(U) .* ϕV[inside]) < abs(U) ? Umax = U : Umin = U
        else
            maximum(ϕρ[inside] .+ T(U) .* ϕV[inside]) - minimum(ϕρ[inside] .+ T(U) .* ϕV[inside]) < abs(U) ? Umin = U : Umax = U
        end        
    end

    if verbose
        @info "The depletion voltage is around $(round(U, digits = Int(ceil(-log10(tolerance))))) ± $(tolerance) V applied to contact $(contact_id)."
        if (potential_range[2] - U) < tolerance || (U - potential_range[1]) < tolerance
            @warn "The depletion voltage is probably not in the specified range $(potential_range.*u"V")."
        end
    end
    return U
end

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