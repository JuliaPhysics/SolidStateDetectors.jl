function _get_abs_bias_voltage(det::SolidStateDetector{T}) where {T <: SSDFloat}
    potentials::Vector{T} = map(c -> c.potential, det.contacts)
    return (maximum(potentials) - minimum(potentials)) * u"V"
end


export calculate_capacitance
@doc raw"""
    calculate_capacitance(sim::Simulation; consider_multiplicity::Bool = true)

Returns the capacitance, ``C``, of a [`SolidStateDetector`](@ref) in a given [`Simulation`](@ref) in units of pF
calculated via 
```math
C = 2 W_{E} / V_{BV}~,
```
where ``W_{E}`` is the energy stored in the electric field ([`calculate_stored_energy`](@ref)) 
created by the detector with the applied bias voltage ``V_{BV}``. 

## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) with `sim.detector` for which the capacitance is calculated.

## Keywords 
* `consider_multiplicity::Bool = true`: Whether symmetries of the system should be taken into account. 
    For example, in case of true coaxial detector center around the origin and calculated on a cartesian grid 
    with the `x-axis` going from `[0, x_max]` and the `y-axis` going from `[0, y_max]` the multiplicity is 4
    and, if `consider_multiplicity == true`, the returned value is already multiplied by 4.

!!! danger 
    In general, this method is only valid if the system, detector and its surroundings, is symmetric with respect to the contacts
    and the impurity density is zero and no fixed charge densities are present (or in the limit of very high bias voltages where
    the contribution of those densities become negligible).
    Then, the returned capacitance equals the capacitance of each contact.
    
    For asymmetric cases and zero impurity/charge densities, 
    the calculated capacitance of this method equals the contact capacitance
    of the contact which is not at ground.

    In general, use [`calculate_capacitance(::Simulation, ::Int)`](@ref) 
    to calculate the capacitance of a specific contact via its weighting potential.
"""
function calculate_capacitance(sim::Simulation; consider_multiplicity::Bool = true)
    W = calculate_stored_energy(sim; consider_multiplicity)
    return uconvert(u"pF", 2 * W / (_get_abs_bias_voltage(sim.detector)^2))
end


@doc raw"""
    calculate_capacitance(sim::Simulation, contact_id::Int = 1; consider_multiplicity::Bool = true)

Returns the capacitance, ``C``, of the contact with ID `contact_id` in units of pF calculated via 
```math
C = 2 W_{WP} / 1 V^2~,
```
where ``W_{WP}`` is the (pseudo) energy stored in the "electric field" (gradient) of the weighting potential of the contact. 

## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) with `sim.detector` for which the capacitance is calculated.
* `contact_id::Int`: The ID of the contact for which the capacitance should be calculated.

## Keywords 
* `consider_multiplicity::Bool = true`: Whether symmetries of the system should be taken into account. 
    For example, in case of true coaxial detector center around the origin and calculated on a cartesian grid 
    with the `x-axis` going from `[0, x_max]` and the `y-axis` going from `[0, y_max]` the multiplicity is 4
    and, if `consider_multiplicity == true`, the returned value is already multiplied by 4.
"""
function calculate_capacitance(sim::Simulation, contact_id::Int; consider_multiplicity::Bool = true)
    @assert !ismissing(sim.weighting_potentials[contact_id]) "Weighting potential of contact $contact_id has not been calculated yet. Please run `calculate_weighting_potential!(sim, $contact_id)` first."
    W = calculate_stored_energy(sim, contact_id; consider_multiplicity)
    return uconvert(u"pF", 2 * W / u"V^2")
end


