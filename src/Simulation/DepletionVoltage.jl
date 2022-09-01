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
    estimate_depletion_voltage( sim::Simulation{T}, contact_id::Int, potential_range::AbstractRange; kwargs... )::T

Estimates the potential (in V) needed at the [`Contact`](@ref) with id `contact_id`
to fully deplete the detector in a given [`Simulation`](@ref). 

## Arguments 
* `sim::Simulation{T}`: [`Simulation`](@ref) for which the depletion voltage should be determined.
* `contact_id::Int`: The `id` of the [`Contact`](@ref) at which the potential is applied.

    
## Keywords
* `field_sim_settings::NamedTuple = (verbose = false,)`: NamedTuple of simulation passed further to the field calculation functions.

## Example 
```julia
using SolidStateDetectors
sim = Simulation(SSD_examples[:InvertedCoax])
SolidStateDetectors.estimate_depletion_voltage(sim, 2, field_sim_settings = (verbose = true,))

```

!!! note
    The accuracy of the result depends on the precision of the initial simulation.
    
See also [`is_depleted`](@ref).
"""
function estimate_depletion_voltage(
        sim::Simulation{T}, 
        contact_id::Int,
        field_sim_settings::NamedTuple = (verbose = false,)
    ) where {T <: AbstractFloat}
    
    simDV = Simulation{T}(sim.config_dict)

    simDV.detector = SolidStateDetector(simDV.detector, contact_id = contact_id, contact_potential = 0)

    max_threads = Base.Threads.nthreads()

    fss = merge((
        convergence_limit = 1e-7,
        max_tick_distance = 2.0u"mm", # this should be like world size / 100 or so
        refinement_limits = [0.2, 0.1, 0.05, 0.025, 0.01], 
        sor_consts = (1.0, 1.0),
        use_nthreads = max_threads > 16 ? 16 : max_threads,
        n_iterations_between_checks = 20, 
        depletion_handling = false
    ), field_sim_settings)

    calculate_electric_potential!(simDV; fss...)
    calculate_weighting_potential!(simDV, contact_id; fss...)
    _adapt_weighting_potential_to_electric_potential_grid!(simDV, contact_id)
    
    grid_size = size(simDV.point_types.grid)
    depletion_voltage_field = zeros(T, grid_size)

    _estimate_depletion_voltage_factor!(
        depletion_voltage_field, 
        grid_size, 
        simDV.point_types, 
        simDV.electric_potential.data, 
        simDV.weighting_potentials[contact_id].data 
    )

    depletion_voltage = maximum(depletion_voltage_field) * u"V"

    # estimated_uncertainty = fss.refinement_limits[end] * depletion_voltage
    return depletion_voltage
end


function _estimate_depletion_voltage_factor(
        center_imp_value::T, 
        center_zero_imp_value::T, 
        neighbor_imp_values::AbstractVector{T}, 
        neighbor_zero_imp_values::AbstractVector{T}
    ) where {T}
    estimation = abs(center_imp_value) / center_zero_imp_value 
    neighbor_values = neighbor_zero_imp_values .* estimation .+ neighbor_imp_values
    vmin, vmax = extrema(neighbor_values)
    center_value = center_zero_imp_value * estimation + center_imp_value
    if center_value < vmin
        estimation += (vmin - center_value) / center_zero_imp_value 
    elseif center_value > vmax
        estimation -= (vmax - center_value) / center_zero_imp_value
    end
    # neighbor_values = neighbor_zero_imp_values .* estimation .+ neighbor_imp_values
    # center_value = center_zero_imp_value * estimation + center_imp_value
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
        grid_size, 
        pts, 
        imp_field, 
        zero_imp_field
    ) where {T}
    @inbounds begin
        for i3 in 1:grid_size[3]
            for i2 in 1:grid_size[2]
                for i1 in 1:grid_size[1]
                    if is_pn_junction_point_type(pts[i1, i2, i3])
                        center_imp_value = imp_field[i1, i2, i3]
                        center_zero_imp_value = zero_imp_field[i1, i2, i3]
                        neighbor_imp_values::SVector{6, T} = _replace_NaN_with_minimum(@SVector T[
                            i1 < grid_size[1] && is_pn_junction_point_type(pts[i1+1, i2, i3]) ? imp_field[i1+1, i2, i3] : T(NaN),
                            i1 > 1 && is_pn_junction_point_type(pts[i1-1, i2, i3]) ? imp_field[i1-1, i2, i3] : T(NaN),
                            i2 < grid_size[2] && is_pn_junction_point_type(pts[i1, i2+1, i3]) ? imp_field[i1, i2+1, i3] : T(NaN),
                            i2 > 1 && is_pn_junction_point_type(pts[i1, i2-1, i3]) ? imp_field[i1, i2-1, i3] : T(NaN),
                            i3 < grid_size[3] && is_pn_junction_point_type(pts[i1, i2, i3+1]) ? imp_field[i1, i2, i3+1] : T(NaN),
                            i3 > 1 && is_pn_junction_point_type(pts[i1, i2, i3-1]) ? imp_field[i1, i2, i3-1] : T(NaN)
                        ])
                        neighbor_zero_imp_values::SVector{6, T} = _replace_NaN_with_minimum(@SVector T[
                            i1 < grid_size[1] && is_pn_junction_point_type(pts[i1+1, i2, i3]) ? zero_imp_field[i1+1, i2, i3] : T(NaN),
                            i1 > 1 && is_pn_junction_point_type(pts[i1-1, i2, i3]) ? zero_imp_field[i1-1, i2, i3] : T(NaN),
                            i2 < grid_size[2] && is_pn_junction_point_type(pts[i1, i2+1, i3]) ? zero_imp_field[i1, i2+1, i3] : T(NaN),
                            i2 > 1 && is_pn_junction_point_type(pts[i1, i2-1, i3]) ? zero_imp_field[i1, i2-1, i3] : T(NaN),
                            i3 < grid_size[3] && is_pn_junction_point_type(pts[i1, i2, i3+1]) ? zero_imp_field[i1, i2, i3+1] : T(NaN),
                            i3 > 1 && is_pn_junction_point_type(pts[i1, i2, i3-1]) ? zero_imp_field[i1, i2, i3-1] : T(NaN)
                        ])

                        depletion_voltage_field[i1, i2, i3] = _estimate_depletion_voltage_factor(
                            center_imp_value, center_zero_imp_value, neighbor_imp_values, neighbor_zero_imp_values
                        )
                    end
                end
            end
        end
    end
end