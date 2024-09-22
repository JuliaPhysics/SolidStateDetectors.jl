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

function Event(location::AbstractCoordinatePoint{T}, energy::RealQuantity = one(T))::Event{T} where {T <: SSDFloat}
    evt = Event{T}()
    evt.locations = VectorOfArrays([[CartesianPoint(location)]])
    evt.energies = VectorOfArrays([[T(to_internal_units(energy))]])
    return evt
end

function Event(locations::Vector{<:AbstractCoordinatePoint{T}}, energies::Vector{<:RealQuantity} = ones(T, length(locations)))::Event{T} where {T <: SSDFloat}
    evt = Event{T}()
    evt.locations = VectorOfArrays(broadcast(pt -> [CartesianPoint(pt)], locations))
    evt.energies = VectorOfArrays(broadcast(E -> [T(to_internal_units(E))], energies))
    return evt
end

function Event(locations::Vector{<:Vector{<:AbstractCoordinatePoint{T}}}, energies::Vector{<:Vector{<:RealQuantity}} = [[one(T) for i in j] for j in locations])::Event{T} where {T <: SSDFloat}
    evt = Event{T}()
    evt.locations = VectorOfArrays(broadcast(pts -> CartesianPoint{T}.(pts), locations))
    evt.energies = VectorOfArrays(Vector{T}.(to_internal_units.(energies)))
    return evt
end

function Event(nbcc::NBodyChargeCloud{T})::Event{T} where {T <: SSDFloat}
    evt = Event{T}()
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
               radius::Vector{<:RealQuantity} = radius_guess.(T.(to_internal_units.(energies)), particle_type)
              )::Event{T} where {T <: SSDFloat, PT <: ParticleType}
    
    return Event(broadcast(i -> 
                NBodyChargeCloud(locations[i], energies[i], N, particle_type, 
                radius = T(to_internal_units(radius[i])), number_of_shells = number_of_shells),
           eachindex(locations)))
end

function Event(evt::NamedTuple{(:evtno, :detno, :thit, :edep, :pos),
         <:Tuple{
            Union{Integer, AbstractVector{<:Integer}},
            Union{Integer, AbstractVector{<:Integer}},
            AbstractVector{<:RealQuantity},
            AbstractVector{<:RealQuantity},
            AbstractVector{<:AbstractVector{<:RealQuantity}}
        }}, T = missing)
    if ismissing(T) T = eltype(to_internal_units(evt[:edep][:])) end

    evt = Event(
        [CartesianPoint{T}.(to_internal_units.(evt[:pos][:]))],
        [T.(to_internal_units.(evt[:edep][:]))]
    )

    return evt
end

in(evt::Event, det::SolidStateDetector) = all( pt -> pt in det, flatview(evt.locations))
in(evt::Event, sim::Simulation) = all( pt -> pt in sim.detector, flatview(evt.locations))

function move_charges_inside_semiconductor!(evt::Event{T}, det::SolidStateDetector{T}; fraction::T = T(0.2))::Event{T} where {T <: SSDFloat}
    move_charges_inside_semiconductor!(evt.locations, evt.energies, det; fraction)
    evt
end

function move_charges_inside_semiconductor!(
        locations::AbstractVector{<:AbstractVector{CartesianPoint{T}}}, energies::AbstractVector{<:AbstractVector{T}},
        det::SolidStateDetector{T}; fraction::T = T(0.2)) where {T <: SSDFloat}
    for n in eachindex(locations)
        idx_in = broadcast( pt -> pt in det.semiconductor, locations[n]);
        if !all(idx_in)
            charge_center = sum(locations[n] .* energies[n]) / sum(energies[n])
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
                                locations[n][m] = (1 - fraction * proj) * i + fraction * proj * charge_center
                            end
                        end
                    end
                end
            end
            charge_center_new = sum(locations[n] .* energies[n]) / sum(energies[n])
            @warn "$(sum(.!idx_in)) charges of the charge cloud at $(round.(charge_center, digits = (T == Float64 ? 12 : 6)))"*
            " are outside. Moving them inside...\nThe new charge center is at $(round.(charge_center_new, digits = (T == Float64 ? 12 : 6))).\n"
        end
    end
    nothing
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
* `verbose = true`: Activate or deactivate additional info output.

## Example 
```julia 
drift_charges!(evt, sim, Δt = 1u"ns", verbose = false)
```

!!! note
    Using values with units for `Δt` requires the package [Unitful.jl](https://github.com/PainterQubits/Unitful.jl).
"""
function drift_charges!(evt::Event{T}, sim::Simulation{T}; max_nsteps::Int = 1000, Δt::RealQuantity = 5u"ns", diffusion::Bool = false, self_repulsion::Bool = false, verbose::Bool = true)::Nothing where {T <: SSDFloat}
    !in(evt, sim.detector) && move_charges_inside_semiconductor!(evt, sim.detector)
    evt.drift_paths = drift_charges(sim, evt.locations, evt.energies, Δt = Δt, max_nsteps = max_nsteps, diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose)
    nothing
end
function get_signal!(evt::Event{T}, sim::Simulation{T}, contact_id::Int; Δt::RealQuantity = 5u"ns")::Nothing where {T <: SSDFloat}
    @assert !ismissing(evt.drift_paths) "No drift path for this event. Use `drift_charges(evt::Event, sim::Simulation)` first."
    @assert !ismissing(sim.weighting_potentials[contact_id]) "No weighting potential for contact $(contact_id). Use `calculate_weighting_potential!(sim::Simulation, contact_id::Int)` first."
    if ismissing(evt.waveforms)
        evt.waveforms = Union{Missing, RadiationDetectorSignals.RDWaveform}[missing for i in eachindex(sim.detector.contacts)]
    end
    evt.waveforms[contact_id] = get_signal(sim, evt.drift_paths, flatview(evt.energies), contact_id, Δt = Δt)
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

## Example 
```julia
get_signals!(evt, sim, Δt = 1u"ns") # if evt.drift_paths were calculated in time steps of 1ns
```

!!! note 
    This method only works if `evt.drift_paths` has already been calculated and is not `missing`.
"""
function get_signals!(evt::Event{T}, sim::Simulation{T}; Δt::RealQuantity = 5u"ns")::Nothing where {T <: SSDFloat}
    @assert !ismissing(evt.drift_paths) "No drift path for this event. Use `drift_charges!(evt::Event, sim::Simulation)` first."
    if ismissing(evt.waveforms)
        evt.waveforms = Union{Missing, RadiationDetectorSignals.RDWaveform}[missing for i in eachindex(sim.detector.contacts)]
    end
    for contact in sim.detector.contacts
        if any(ismissing, sim.weighting_potentials) "No weighting potential(s) for some contact(s).." end
        if !ismissing(sim.weighting_potentials[contact.id])
            evt.waveforms[contact.id] = get_signal(sim, evt.drift_paths, flatview(evt.energies), contact.id, Δt = Δt)
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
* `verbose = true`: Activate or deactivate additional info output.

## Example 
```julia
simulate!(evt, sim, Δt = 1u"ns", verbose = false)
```

See also [`drift_charges!`](@ref) and [`get_signals!`](@ref).
"""
function simulate!(evt::Event{T}, sim::Simulation{T}; max_nsteps::Int = 1000, Δt::RealQuantity = 5u"ns", diffusion::Bool = false, self_repulsion::Bool = false, verbose::Bool = true)::Nothing where {T <: SSDFloat}
    drift_charges!(evt, sim, max_nsteps = max_nsteps, Δt = Δt, diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose)
    get_signals!(evt, sim, Δt = Δt)
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
    new_times = if TV <: AbstractRange
        range( zero(first(wv.time)), step = step(wv.time), length = total_waveform_length )
    else
        error("Not yet definted for timestamps of type `$(TV)`")
    end
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
function get_electron_and_hole_contribution(evt::Event{T}, sim::Simulation{T, S}, contact_id::Int
            )::NamedTuple{(:electron_contribution, :hole_contribution), <:Tuple{RDWaveform, RDWaveform}} where {T <: SSDFloat, S}
    
    @assert !ismissing(evt.drift_paths) "The charge drift is not yet simulated. Please use `drift_charges!(evt, sim)`!"
    
    dt::T = T(ustrip(u"ns", diff(evt.drift_paths[1].timestamps_e)[1]*u"s"))
    wp::Interpolations.Extrapolation{T, 3} = interpolated_scalarfield(sim.weighting_potentials[contact_id])
    signal_e::Vector{T} = zeros(T, length(maximum(map(p -> p.timestamps_e, evt.drift_paths))))
    signal_h::Vector{T} = zeros(T, length(maximum(map(p -> p.timestamps_h, evt.drift_paths))))

    ctm = sim.detector.semiconductor.charge_trapping_model
    for i in eachindex(evt.drift_paths)
        energy = flatview(evt.energies)[i]
        
        dp_e::Vector{CartesianPoint{T}} = evt.drift_paths[i].e_path
        dp_e_t::Vector{T} = evt.drift_paths[i].timestamps_e
        add_signal!(signal_e, dp_e_t, dp_e, dp_e_t, -energy, wp, S, ctm)
        
        dp_h::Vector{CartesianPoint{T}} = evt.drift_paths[i].h_path
        dp_h_t::Vector{T} = evt.drift_paths[i].timestamps_h
        add_signal!(signal_h, dp_h_t, dp_h, dp_h_t, energy, wp, S, ctm)
    end
    unitless_energy_to_charge = _convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)
    return (electron_contribution = RDWaveform(range(zero(T) * u"ns", step = dt * u"ns", length = length(signal_e)), signal_e * unitless_energy_to_charge),
            hole_contribution = RDWaveform(range(zero(T) * u"ns", step = dt * u"ns", length = length(signal_h)), signal_h * unitless_energy_to_charge))
end

export get_electron_and_hole_contribution


"""
    plot_electron_and_hole_contribution(evt::Event{T}, sim::Simulation{T}, contact_id::Int)
    
Plots the waveform as well as the electron and hole contribution to the waveform 
of a [`Contact`](@ref) with a given `contact_id` of an [`Event`](@ref).

## Arguments
* `evt::Event{T}`: [`Event`](@ref) in which the charges have already been drifted.
* `sim::Simulation{T}`: [`Simulation`](@ref) which defines the setup in which the charges in `evt` were drifted.
* `contact_id::Int`: The `id` of the [`Contact`](@ref) for which the waveform should be split into electron and hole contribution.

## Keywords 
* `n_samples::Int`: Number of samples with which the waveforms will be plotted. The default is the number of samples of the original waveform.

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
plot_electron_and_hole_contribution(evt, sim, contact_id, n_samples = 300)
```

!!! note 
    The drift paths in `evt` need to be calculated using [`drift_charges!`](@ref) before calling this function.
    
!!! note 
    This method requires to load the package [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
    
See also [`get_electron_and_hole_contribution`](@ref).
"""
function plot_electron_and_hole_contribution end

@userplot Plot_electron_and_hole_contribution
@recipe function f(gdd::Plot_electron_and_hole_contribution; linewidth = 2, n_samples = missing)
    
    sim::Simulation = gdd.args[2]
    T = get_precision_type(sim)
    evt::Event{T} = gdd.args[1]
    contact_id::Int = gdd.args[3]
    
    ismissing(n_samples) ? n_samples = length(evt.waveforms[contact_id].signal) : nothing
    
    wf::NamedTuple{(:electron_contribution, :hole_contribution), <:Tuple{RDWaveform, RDWaveform}} = get_electron_and_hole_contribution(evt, sim, contact_id)
    
    unitformat --> :slash
    yguide --> "Charge"

    @series begin
        linecolor := :red 
        linewidth --> linewidth
        label --> "Electron contribution"
        add_baseline_and_extend_tail(wf.electron_contribution, 0, n_samples)
    end
    
    @series begin
        linecolor := :green 
        linewidth --> linewidth
        label --> "Hole contribution"
        add_baseline_and_extend_tail(wf.hole_contribution, 0, n_samples)
    end
    
    @series begin
        linecolor --> contact_id
        linewidth --> linewidth
        seriesalpha --> 0.5
        label --> "Sum"
        add_baseline_and_extend_tail(evt.waveforms[contact_id], 0, n_samples)
    end
end