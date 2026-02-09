"""
  mutable struct Event{T <: SSDFloat}
         
Collection struct for individual events. 
This (mutable) struct is meant to be used to look at individual events,
not to process a huge amount of events.

## Fields
* `locations::VectorOfArrays{CartesianPoint{T}}`: Vector of the positions of all hits of the event.
* `energies::VectorOfArrays{T}`: Vector of energies corresponding to the hits of the event.
* `drift_paths::Union{Vector{EHDriftPath{T}}, Missing}`: Calculated drift paths of each hit position. 
* `waveforms::Union{Vector{<:Any}, Missing}`: Generated signals (waveforms) of the event.
"""
mutable struct Event{T <: SSDFloat}
    locations::VectorOfArrays{CartesianPoint{T}}
    energies::VectorOfArrays{T}
    drift_paths::Union{Vector{EHDriftPath{T}}, Missing}
    waveforms::Union{Vector{<:Any}, Missing}
    Event{T}() where {T <: SSDFloat} = new{T}(VectorOfArrays(Vector{CartesianPoint{T}}[]), VectorOfArrays(Vector{T}[]), missing, missing)
end

function Event(location::AbstractCoordinatePoint, energy::RealQuantity = 1)
    pt_internal = to_internal_units(location)
    T = eltype(pt_internal)
    evt = Event{T}()
    evt.locations = VectorOfArrays([[pt_internal]])
    evt.energies  = VectorOfArrays([[T(to_internal_units(energy))]])
    return evt
end

function Event(locations::Vector{<:AbstractCoordinatePoint}, energies::Vector{<:RealQuantity} = fill(1, length(locations)); max_interaction_distance::Union{<:Real, <:LengthQuantity} = NaN)
    pts_internal = to_internal_units.(locations)
    T = eltype(first(pts_internal))
    d::T = T(to_internal_units(max_interaction_distance))
    @assert isnan(d) || d >= 0 "Max. interaction distance must be positive or NaN (no grouping), but $(max_interaction_distance) was given."
    evt = Event{T}()
    
    if isnan(d) # default: no grouping, the charges from different hits drift independently
        evt.locations = VectorOfArrays(broadcast(pt -> [pt], pts_internal))
        evt.energies  = VectorOfArrays(broadcast(E -> [T(to_internal_units(E))], energies))
    else
        evt.locations, evt.energies = group_points_by_distance(pts_internal, T.(to_internal_units.(energies)), d)
    end
    return evt
end

function Event(locations::Vector{<:Vector{<:AbstractCoordinatePoint}}, energies::Vector{<:Vector{<:RealQuantity}} = [[1 for _ in group] for group in locations])

    T = eltype(to_internal_units(first(first(locations))))
    locs = [[convert(CartesianPoint{T}, to_internal_units(pt)) for pt in group] for group in locations]
    ens  = [[T(to_internal_units(E)) for E in group] for group in energies]
    
    evt = Event{T}()
    evt.locations = VectorOfArrays(locs)
    evt.energies  = VectorOfArrays(ens)
    return evt
end

function Event(nbcc::NBodyChargeCloud{T})::Event{T} where {T <: SSDFloat}
    evt = Event{T}()
    # charges in one NBodyChargeCloud should see each other by default
    evt.locations = VectorOfArrays([nbcc.locations])
    evt.energies = VectorOfArrays([nbcc.energies])
    return evt
end

function Event(nbccs::Vector{<:NBodyChargeCloud{T}})::Event{T} where {T <: SSDFloat}
    evt = Event{T}()
    evt.locations = VectorOfArrays(broadcast(nbcc -> nbcc.locations, nbccs))
    evt.energies = VectorOfArrays(broadcast(nbcc -> nbcc.energies, nbccs))
    return evt
end

function Event(locations::Vector{<:AbstractCoordinatePoint{T}}, energies::Vector{<:RealQuantity}, N::Int; 
        particle_type::Type{PT} = Gamma, number_of_shells::Int = 2,
        radius::Vector{<:Union{<:Real, <:LengthQuantity}} = radius_guess.(T.(to_internal_units.(energies)), particle_type),
        max_interaction_distance::Union{<:Real, <:LengthQuantity} = NaN
    )::Event{T} where {T <: SSDFloat, PT <: ParticleType}

    @assert eachindex(locations) == eachindex(energies) == eachindex(radius)

    d::T = T(to_internal_units(max_interaction_distance))
    @assert isnan(d) || d >= 0 "Max. interaction distance must be positive or NaN (no grouping), but $(max_interaction_distance) was given."

    if isnan(d) # default: no grouping, the charges from different hits drift independently
        Event(
            broadcast(i -> 
                NBodyChargeCloud(locations[i], energies[i], N, 
                radius = T(to_internal_units(radius[i])), number_of_shells = number_of_shells),
            eachindex(locations))
        )
    else # otherwise: group the CENTERS by distance and compose the NBodyChargeClouds around them
        loc, ene, rad = group_points_by_distance(CartesianPoint.(locations), T.(to_internal_units.(energies)), T.(to_internal_units.(radius)), d)
        nbccs = broadcast(i -> 
            broadcast(j -> 
                let nbcc = NBodyChargeCloud(loc[i][j], ene[i][j], N, 
                    radius = rad[i][j], number_of_shells = number_of_shells)
                    nbcc.locations, nbcc.energies
                end, 
            eachindex(loc[i])), 
        eachindex(loc))

        Event(
            broadcast(nbcc -> vcat(getindex.(nbcc, 1)...), nbccs),
            broadcast(nbcc -> vcat(getindex.(nbcc, 2)...), nbccs)
        )
    end
end

function Event(
    locations::Vector{<:Vector{<:AbstractCoordinatePoint}},
    energies::Vector{<:Vector{<:RealQuantity}}, N::Int;
    particle_type::Type = Gamma, number_of_shells::Int = 2,
    radius::Vector{<:Vector{<:RealQuantity}} = map(e -> radius_guess.(to_internal_units.(e), particle_type), energies)
)
    locs_internal = [to_internal_units.(loc) for loc in locations]

    @assert eachindex(locs_internal) == eachindex(energies) == eachindex(radius)
    events = map(i -> Event(locs_internal[i], energies[i], N; particle_type, number_of_shells, radius = radius[i]), eachindex(locations))

    Event(flatview.(getfield.(events, :locations)), flatview.(getfield.(events, :energies)))
end

function Event(evt::NamedTuple{DHE_column_names, <:Tuple{DHE_column_types...}}, T = eltype(to_internal_units(evt.edep[:])))::Event{T}
    return Event(
        [CartesianPoint{T}.(to_internal_units.(evt.pos[:]))],
        [T.(to_internal_units.(evt.edep[:]))];
    )
end

in(evt::Event, det::SolidStateDetector) = all( pt -> pt in det, flatview(evt.locations))
in(evt::Event, sim::Simulation) = in(evt, sim.detector)

function move_charges_inside_semiconductor!(evt::Event{T}, det::SolidStateDetector{T}; fraction::T = T(0.2), verbose::Bool = true) where {T <: SSDFloat}
    move_charges_inside_semiconductor!(evt.locations, evt.energies, det; fraction, verbose)
end

function move_charges_inside_semiconductor!(
        locations::AbstractVector{<:AbstractVector{CartesianPoint{T}}}, energies::AbstractVector{<:AbstractVector{T}},
        det::SolidStateDetector{T}; fraction::T = T(0.2), verbose::Bool = true) where {T <: SSDFloat}
    for n in eachindex(locations)
        idx_in = broadcast( pt -> pt in det.semiconductor, locations[n]);
        if !all(idx_in)
            edep_weights = ustrip.(energies[n])
            charge_center = barycenter(locations[n], Weights(edep_weights))
            @assert charge_center in det.semiconductor "The center of the charge cloud ($(charge_center)) is not inside the semiconductor."
            surf = ConstructiveSolidGeometry.surfaces(det.semiconductor.geometry)
            for (k,m) in enumerate(findall(.!idx_in))
                l = locations[n][m]
                line = ConstructiveSolidGeometry.Line(charge_center, CartesianVector{T}(l - charge_center))
                for s in surf
                    int = ConstructiveSolidGeometry.intersection(s, line)
                    for i in int
                        proj::T = (i - charge_center) ⋅ (l - charge_center) / sum((l - charge_center).^2)
                        if 0 ≤ proj ≤ 1
                            if i in det.semiconductor
                                locations[n][m] = cartesian_zero + (1 - fraction * proj) * (i - cartesian_zero) + fraction * proj * (charge_center - cartesian_zero)
                            end
                        end
                    end
                end
            end
            charge_center_new = barycenter(locations[n], Weights(edep_weights))
            if verbose
                @warn "$(sum(.!idx_in)) charges of the charge cloud at $(geom_round(charge_center))"*
                " are outside. Moving them inside...\nThe new charge center is at $(geom_round(charge_center_new)).\n"
            end
        end
    end
    return nothing
end

"""
    drift_charges!(evt::Event{T}, sim::Simulation{T}; kwargs...)::Nothing where {T <: SSDFloat}

Calculates the electron and hole drift paths for the given [`Event`](@ref) and [`Simulation`](@ref)
    and stores them in `evt.drift_paths`.
    
## Arguments 
* `evt::Event{T}`: [`Event`](@ref) for which the charges should be drifted.
* `sim::Simulation{T}`: [`Simulation`](@ref) which defines the setup in which the charges in `evt` should drift.

## Keywords
* `max_nsteps::Int = 1000`: Maximum number of steps in the drift of each hit. 
* `Δt::RealQuantity = 5u"ns"`: Time step used for the drift.
* `diffusion::Bool = false`: Activate or deactive diffusion of charge carriers via random walk.
* `self_repulsion::Bool = false`: Activate or deactive self-repulsion of charge carriers of the same type.
* `end_drift_when_no_field::Bool = true`: Activate or deactive drifting termination when the electric field is exactly zero.
* `geometry_check::Bool = false`: Perform extra geometry checks when determining if charge carriers have reached a contact.
* `verbose = true`: Activate or deactivate additional info output.

## Example 
```julia 
drift_charges!(evt, sim, Δt = 1u"ns", verbose = false)
```

!!! note
    Using values with units for `Δt` requires the package [Unitful.jl](https://github.com/JuliaPhysics/Unitful.jl).
"""
function drift_charges!(evt::Event{T}, sim::Simulation{T}; max_nsteps::Int = 1000, Δt::RealQuantity = 5u"ns", diffusion::Bool = false, self_repulsion::Bool = false, end_drift_when_no_field::Bool = true, geometry_check::Bool = false, verbose::Bool = true)::Nothing where {T <: SSDFloat}
    !in(evt, sim) && move_charges_inside_semiconductor!(evt, sim.detector; verbose)
    evt.drift_paths = drift_charges(sim, evt.locations, evt.energies; Δt, max_nsteps, diffusion, self_repulsion, end_drift_when_no_field, geometry_check, verbose)
    nothing
end
function get_signal!(evt::Event{T}, sim::Simulation{T}, contact_id::Int; Δt::RealQuantity = 5u"ns", signal_unit::Unitful.Units = u"e_au")::Nothing where {T <: SSDFloat}
    ismissing(evt.drift_paths) && throw(AssertionError("The charge drift is not yet simulated. Please run `drift_charges!(evt::Event, sim::Simulation)` first."))
    ismissing(sim.weighting_potentials[contact_id]) && throw(AssertionError("No weighting potential for contact $(contact_id). Please run `calculate_weighting_potential!(sim::Simulation, contact_id::Int)` first."))
    if ismissing(evt.waveforms)
        evt.waveforms = Union{Missing, RadiationDetectorSignals.RDWaveform}[missing for i in eachindex(sim.detector.contacts)]
    end
    evt.waveforms[contact_id] = get_signal(sim, evt.drift_paths, flatview(evt.energies), contact_id; Δt, signal_unit)
    nothing
end


"""
    get_signals!(evt::Event{T}, sim::Simulation{T}; Δt::RealQuantity = 5u"ns")::Nothing where {T <: SSDFloat}
    
Generates the signals/waveforms from the drift paths of an [`Event`](@ref) for each [`Contact`](@ref),
for which a [`WeightingPotential`](@ref) is specified in `sim.weighting_potentials`.
    
The output is stored in `evt.waveforms`.

## Arguments 
* `evt::Event{T}`: [`Event`](@ref) for which the waveforms should be generated.
* `sim::Simulation{T}`: [`Simulation`](@ref) which defines the setup in which the waveforms are generated.

## Keywords
* `Δt::RealQuantity = 5u"ns"`: Time steps with which the drift paths were calculated.
* `signal_unit::Unitful.Units = u"e_au"`: Unit of the returned waveform (charge or energy).

## Example 
```julia
SolidStateDetectors.get_signals!(evt, sim, Δt = 1u"ns") # if evt.drift_paths were calculated in time steps of 1ns
```

!!! note 
    This method only works if `evt.drift_paths` has already been calculated and is not `missing`.
"""
function get_signals!(evt::Event{T}, sim::Simulation{T}; kwargs...)::Nothing where {T <: SSDFloat}
    any(ismissing, sim.weighting_potentials) && @warn "No weighting potential(s) for some contact(s).."
    for contact in sim.detector.contacts
        if !ismissing(sim.weighting_potentials[contact.id])
            get_signal!(evt, sim, contact.id; kwargs...)
        end
    end
    nothing
end

"""
    simulate!(evt::Event{T}, sim::Simulation{T}; kwargs...)::Nothing where {T <: SSDFloat}

Simulates the waveforms for the [`Event`](@ref) for a given [`Simulation`](@ref) by

1. calculating the drift paths of all energy hits, at `evt.locations` and 
2. generating the waveforms for each [`Contact`](@ref), for which a [`WeightingPotential`](@ref) is specified in `sim.weighting_potentials`.

The output is stored in `evt.drift_paths` and `evt.waveforms`.

## Arguments 
* `evt::Event{T}`: [`Event`](@ref) for which the charges should be drifted.
* `sim::Simulation{T}`: [`Simulation`](@ref) which defines the setup in which the charges in `evt` should drift.

## Keywords
* `max_nsteps::Int = 1000`: Maximum number of steps in the drift of each hit. 
* `Δt::RealQuantity = 5u"ns"`: Time step used for the drift.
* `diffusion::Bool = false`: Activate or deactive diffusion of charge carriers via random walk.
* `self_repulsion::Bool = false`: Activate or deactive self-repulsion of charge carriers of the same type.
* `signal_unit::Unitful.Units = u"e_au"`: Unit of the returned waveform (charge or energy).
* `end_drift_when_no_field::Bool = true`: Activate or deactive drifting termination when the electric field is exactly zero.
* `geometry_check::Bool = false`: Perform extra geometry checks when determining if charge carriers have reached a contact.
* `verbose = true`: Activate or deactivate additional info output.

## Example 
```julia
simulate!(evt, sim, Δt = 1u"ns", verbose = false)
```

See also [`drift_charges!`](@ref) and [`get_signals!`](@ref).
"""
function simulate!(evt::Event{T}, sim::Simulation{T}; max_nsteps::Int = 1000, Δt::RealQuantity = 5u"ns", diffusion::Bool = false, self_repulsion::Bool = false, signal_unit::Unitful.Units = u"e_au", end_drift_when_no_field::Bool = true, geometry_check::Bool = false, verbose::Bool = true)::Nothing where {T <: SSDFloat}
    drift_charges!(evt, sim; max_nsteps, Δt, diffusion, self_repulsion, end_drift_when_no_field, geometry_check, verbose)
    get_signals!(evt, sim; Δt, signal_unit)
    nothing
end


"""
    add_baseline_and_extend_tail(wv::RadiationDetectorSignals.RDWaveform{T,U,TV,UV}, n_baseline_samples::Int, total_waveform_length::Int) where {T,U,TV,UV}

Adds a zero-valued baseline in front of the waveform `wv` and extends (or cuts off) the waveform at the end with the last value of `wv`.
A waveform of length `total_waveform_length` is returned.

## Arguments
* `wv::RadiationDetectorSignals.RDWaveform{T,U,TV,UV}`: A waveform (signal over time).
* `n_baseline_samples::Int`: Number of samples added in front of the waveform with values 0. 
* `total_waveform_length::Int`: Number of samples of the extended waveform which is returned.

## Examples
    add_baseline_and_extend_tail(wv, 1000, 5000)   

!!! note 
    This functions assumes that the time steps between the samples of the input waveform `wv` are the same. Thus, that the input waveform is
    sampled with a fixed frequency. 
"""
function add_baseline_and_extend_tail(wv::RadiationDetectorSignals.RDWaveform{T,U,TV,UV}, n_baseline_samples::Int, total_waveform_length::Int) where {T,U,TV,UV}
    new_signal::Vector{eltype(UV)} = Vector{eltype(UV)}(undef, total_waveform_length)
    new_signal[1:n_baseline_samples] .= zero(eltype(UV))
    if length(wv.signal) <= total_waveform_length - n_baseline_samples
        new_signal[n_baseline_samples+1:n_baseline_samples+length(wv.signal)] = wv.signal
        new_signal[n_baseline_samples+length(wv.signal)+1:end] .= wv.signal[end]
    else
        new_signal[n_baseline_samples+1:end] = wv.signal[1:total_waveform_length - n_baseline_samples]
    end
    TV <: AbstractRange || throw(ArgumentError("Not yet defined for timestamps of type `$(TV)`"))
    new_times = range(first(wv.time), step = step(wv.time), length = total_waveform_length )
    return RDWaveform( new_times, new_signal )
end


"""
    get_electron_and_hole_contribution(evt::Event{T}, sim::Simulation{T}, contact_id::Int)
    
Returns the electron and hole contribution to the waveform of a [`Contact`](@ref) with a given
`contact_id` of an [`Event`](@ref) as a `NamedTuple` with two entries: 
`electron_contribution` and `hole_contribution`.

## Arguments
* `evt::Event{T}`: [`Event`](@ref) in which the charges have already been drifted.
* `sim::Simulation{T}`: [`Simulation`](@ref) which defines the setup in which the charges in `evt` were drifted.
* `contact_id::Int`: The `id` of the [`Contact`](@ref) for which the waveform should be split into electron and hole contribution.

## Example
```julia 
using Plots
using SolidStateDetector
T = Float32

simulation = Simulation{T}(SSD_examples[:InvertedCoax])
simulate!(simulation)
event = Event([CartesianPoint{T}(0.02,0.01,0.05)])
simulate!(event, simulation)

contact_id = 1
wf = get_electron_and_hole_contribution(evt, sim, contact_id)
```

!!! note 
    The drift paths in `evt` need to be calculated using [`drift_charges!`](@ref) before calling this function.
    
See also [`plot_electron_and_hole_contribution`](@ref).
"""
function get_electron_and_hole_contribution(evt::Event{T}, sim::Simulation{T, S}, contact_id::Int; signal_unit::Unitful.Units = u"e_au",
            )::NamedTuple{(:electron_contribution, :hole_contribution), <:Tuple{RDWaveform, RDWaveform}} where {T <: SSDFloat, S}
    
    ismissing(evt.drift_paths) && throw(AssertionError("The charge drift is not yet simulated. Please run `drift_charges!(evt::Event, sim::Simulation)` first."))
    
    dt::T = T(ustrip(u"ns", diff(evt.drift_paths[1].timestamps_e)[1]*u"s"))
    wp::Interpolations.Extrapolation{T, 3} = interpolated_scalarfield(sim.weighting_potentials[contact_id])
    signal_e::Vector{T} = zeros(T, length(maximum(map(p -> p.timestamps_e, evt.drift_paths))))
    signal_h::Vector{T} = zeros(T, length(maximum(map(p -> p.timestamps_h, evt.drift_paths))))

    ctm = sim.detector.semiconductor.charge_trapping_model
    for i in eachindex(evt.drift_paths)
        energy = flatview(evt.energies)[i]
        
        dp_e::Vector{CartesianPoint{T}} = evt.drift_paths[i].e_path
        dp_e_t::Vector{T} = evt.drift_paths[i].timestamps_e
        add_signal!(signal_e, dp_e_t, dp_e, dp_e_t, -energy, wp, sim.point_types, ctm)
        
        dp_h::Vector{CartesianPoint{T}} = evt.drift_paths[i].h_path
        dp_h_t::Vector{T} = evt.drift_paths[i].timestamps_h
        add_signal!(signal_h, dp_h_t, dp_h, dp_h_t, energy, wp, sim.point_types, ctm)
    end
    calibration_factor::Quantity{T, dimension(signal_unit)} = _convert_internal_energy_to_external_unit(signal_unit, sim.detector.semiconductor.material)
    return (electron_contribution = RDWaveform(range(zero(T) * u"ns", step = dt * u"ns", length = length(signal_e)), signal_e * calibration_factor),
            hole_contribution = RDWaveform(range(zero(T) * u"ns", step = dt * u"ns", length = length(signal_h)), signal_h * calibration_factor))
end

export get_electron_and_hole_contribution
