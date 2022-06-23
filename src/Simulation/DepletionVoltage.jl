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
    get_depletion_voltage( sim::Simulation{T}, contact_id::Int, potential_range::AbstractRange; kwargs... )::T

Returns the potential (in V) needed at the [`Contact`](@ref) with id `contact_id`
to fully deplete the detector in a given [`Simulation`](@ref). For this, all other
contact potentials are set to `0` and the potential at the specified contact is
increased or decreased according to the `potential_range`. The depletion voltage
is set to the potential for which a previously undepleted detector becomes depleted,
or for which a previously depleted detector becomes undepleted.

## Arguments 
* `sim::Simulation{T}`: [`Simulation`](@ref) for which the depletion voltage should be determined.
* `contact_id::Int`: The `id` of the [`Contact`](@ref) at which the potential is applied.
* `potential_range::AbstractRange`: Range of potentials to be tested.
    
## Keywords
* `verbose::Bool = true`: Activate or deactivate additional info output. Default is `true`.

## Example 
```julia
using SolidStateDetectors
sim = Simulation(SSD_examples[:InvertedCoax])
calculate_electric_potential!(sim)
get_depletion_voltage(sim, 2, 1600:1:2500)
```

!!! warn
    This method only works if the initial `sim` was calculated for a case in which the detector is fully depleted.

!!! note
    The accuracy of the result depends on the precision of the initial simulation.
    
See also [`is_depleted`](@ref).
"""
function get_depletion_voltage(sim::Simulation{T}, contact_id::Int,
            potential_range::AbstractRange = range(extrema(broadcast(c -> c.potential, sim.detector.contacts))..., length = 1001);
            verbose = true)::T where {T <: AbstractFloat}
    
    @assert is_depleted(sim.point_types) "This method only works for fully depleted simulations. Please increase the potentials in the configuration file to a greater value."
    
    if verbose
        @info "Looking for the depletion voltage applied to contact $(contact_id) "*
              "in the range $(extrema(potential_range).*u"V") in steps of $(step(potential_range)*u"V")."
    end
    
    for c in sim.detector.contacts
        if c.potential == 0 && c.id != contact_id continue end
        if ismissing(sim.weighting_potentials[c.id]) || sim.weighting_potentials[c.id].grid != sim.electric_potential.grid
            @info "The weighting potential $(c.id) need to be defined on the same grid as the electric potential. Adjusting weighting potential $(c.id) now."
            _adapt_weighting_potential_to_electric_potential_grid!(sim, c.id)
        end
    end
    
    # ϕρ is the electric potential resulting only from the impurity density
    # and setting all potentials on the contacts to zero
    ϕρ = deepcopy(sim.electric_potential.data)
    for c in sim.detector.contacts
        if c.potential == 0 continue end
        ϕρ .-= c.potential * sim.weighting_potentials[c.id].data
    end
    
    inside = findall(p -> p & pn_junction_bit > 0 && p & update_bit > 0, sim.point_types.data)
    
    ϕmin, ϕmax = extrema((ϕρ .+ potential_range[1] * sim.weighting_potentials[contact_id].data)[inside])
    initial_depletion::Bool = ϕmax - ϕmin < abs(potential_range[1])
    depletion_voltage::T = NaN
    
    I1 = oneunit(CartesianIndex(0,0,0))
    Ifirst, Ilast = first(inside), last(inside)
    scale = zeros(size(sim.electric_potential.data))
    for idx in inside
        for scale_loc in potential_range
            center_pot = T(scale_loc) * sim.weighting_potentials[contact_id].data[idx] + ϕρ[idx]
            min_pot = sim.weighting_potentials[contact_id].data[idx+CartesianIndex(1,0,0)] + ϕρ[idx+CartesianIndex(1,0,0)]
            max_pot = min_pot
            for n_idx in (max(Ifirst,(idx-I1)):min(Ilast,(idx-I1)))
                if n_idx == idx continue end
                local_pot = T(scale_loc) * sim.weighting_potentials[contact_id].data[n_idx] + ϕρ[n_idx]
                if local_pot <=min_pot
                    min_pot = local_pot
                elseif local_pot>=max_pot
                    max_pot = local_pot
                end
            end
            if min_pot<=center_pot<=max_pot  
                scale[idx] = T(scale_loc)
                println(scale_loc)
                break
            end         
        end
    end
    print(scale[inside])
    depletion_voltage = maximum(scale[inside])
    # @showprogress for U in potential_range
    #     ϕmin, ϕmax = extrema((ϕρ .+ T(U) * sim.weighting_potentials[contact_id].data)[inside])
    #     depleted = ϕmax - ϕmin < abs(U)
    #     if (initial_depletion && !depleted) || (!initial_depletion && depleted)
    #         depletion_voltage = T(U)
    #         break
    #     end
    # end
    
    if verbose
        if !isnan(depletion_voltage)
            @info "The depletion voltage of the detector is ($(depletion_voltage) ± $(T(step(potential_range)))) V applied to contact $(contact_id)."
        else
            @warn "The depletion voltage is not in the specified range $(extrema(potential_range).*u"V")."
        end
    end
    return depletion_voltage
end