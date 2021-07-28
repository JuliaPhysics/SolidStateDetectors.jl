"""
    mutable struct Event{T <: SSDFloat}

Collection struct for individual events. 
This (mutable) struct is meant to be used to look at individual events,
not to process a huge amount of events.

## Fields
* `locations::Union{Vector{<:AbstractCoordinatePoint{T}}, Missing}`: Vector of the positions of all hits of the event.
* `energies::Union{Vector{T}, Missing}`: Vector of energies corresponding to the hits of the event.
* `drift_paths::Union{Vector{EHDriftPath{T}}, Missing}`: Calculated drift paths of each hit position. 
* `waveforms::Union{Vector{<:Any}, Missing}`: Generated signals (waveforms) of the event.
"""
mutable struct Event{T <: SSDFloat}
    locations::Union{Vector{<:AbstractCoordinatePoint{T}}, Missing}
    energies::Union{Vector{T}, Missing}
    drift_paths::Union{Vector{EHDriftPath{T}}, Missing}
    waveforms::Union{Vector{<:Any}, Missing}
    Event{T}() where {T <: SSDFloat} = new{T}(missing, missing, missing, missing)
end

function Event(locations::Vector{<:AbstractCoordinatePoint{T}}, energies::Vector{<:RealQuantity{T}} = ones(T, length(locations)))::Event{T} where {T <: SSDFloat}
    evt = Event{T}()
    evt.locations = locations
    evt.energies = to_internal_units(energies)
    evt.waveforms = missing
    return evt
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
        CartesianPoint{T}.(to_internal_units.(evt[:pos][:])),
        T.(to_internal_units.(evt[:edep][:]))
    )

    return evt
end

in(evt::Event, det::SolidStateDetector) = all( pt -> pt in det, evt.locations)
in(evt::Event, sim::Simulation) = all( pt -> pt in sim.detector, evt.locations)

"""
    drift_charges!(evt::Event{T}, sim::Simulation{T}; max_nsteps::Int = 1000, Δt::RealQuantity = 5u"ns", verbose::Bool = true)::Nothing where {T <: SSDFloat}

Calculate the drift paths for the given `evt` and Simulation `sim`
    and stores them under `evt.drift_paths`. 

# Keywords
- `Δt::RealQuantity = 5u"ns"`: Time difference between two time stamps of the drift.
- `max_nsteps::Int = 1000`: Maximum number of steps in the drift of each hit. 
- `verbose = false`: Activate or deactivate additional info output. 
"""
function drift_charges!(evt::Event{T}, sim::Simulation{T}; max_nsteps::Int = 1000, Δt::RealQuantity = 5u"ns", verbose::Bool = true)::Nothing where {T <: SSDFloat}
    evt.drift_paths = drift_charges(sim, CartesianPoint.(evt.locations), Δt = Δt, max_nsteps = max_nsteps, verbose = verbose)
    nothing
end
function get_signal!(evt::Event{T}, sim::Simulation{T}, contact_id::Int; Δt::RealQuantity = 5u"ns")::Nothing where {T <: SSDFloat}
    @assert !ismissing(evt.drift_paths) "No drift path for this event. Use `drift_charges(evt::Event, sim::Simulation)` first."
    @assert !ismissing(sim.weighting_potentials[contact_id]) "No weighting potential for contact $(contact_id). Use `calculate_weighting_potential!(sim::Simulation, contact_id::Int)` first."
    if ismissing(evt.waveforms)
        evt.waveforms = Union{Missing, RadiationDetectorSignals.RDWaveform}[missing for i in eachindex(sim.detector.contacts)]
    end
    evt.waveforms[contact_id] = get_signal(sim, evt.drift_paths, evt.energies, contact_id, Δt = Δt)
    nothing
end
function get_signals!(evt::Event{T}, sim::Simulation{T}; Δt::RealQuantity = 5u"ns")::Nothing where {T <: SSDFloat}
    @assert !ismissing(evt.drift_paths) "No drift path for this event. Use `drift_charges(evt::Event, sim::Simulation)` first."
    if ismissing(evt.waveforms)
        evt.waveforms = Union{Missing, RadiationDetectorSignals.RDWaveform}[missing for i in eachindex(sim.detector.contacts)]
    end
    for contact in sim.detector.contacts
        if any(ismissing, sim.weighting_potentials) "No weighting potential(s) for some contact(s).." end
        if !ismissing(sim.weighting_potentials[contact.id])
            evt.waveforms[contact.id] = get_signal(sim, evt.drift_paths, evt.energies, contact.id, Δt = Δt)
        end
    end
    nothing
end

"""
    simulate!(evt::Event{T}, sim::Simulation{T}; max_nsteps::Int = 1000, Δt::RealQuantity = 5u"ns", verbose::Bool = true)::Nothing where {T <: SSDFloat}

Simulates the `evt` for a given Simulation `sim`:
    1. Simulate the drift paths of all hits, at `evt.locations`, and stores them under `evt.drift_paths`. 
    2. Generate the signal (waveform) for all channels for which a weighting potential is specified in `sim.weighting_potentials`
        and stores them under `evt.waveforms`.

# Keywords
- `Δt::RealQuantity = 5u"ns"`: Time difference between two time stamps of the drift and the signals.
- `max_nsteps::Int = 1000`: Maximum number of steps in the drift of each hit. 
- `verbose = false`: Activate or deactivate additional info output. 
"""
function simulate!(evt::Event{T}, sim::Simulation{T}; max_nsteps::Int = 1000, Δt::RealQuantity = 5u"ns", verbose::Bool = true)::Nothing where {T <: SSDFloat}
    drift_charges!(evt, sim, max_nsteps = max_nsteps, Δt = Δt, verbose = verbose)
    get_signals!(evt, sim, Δt = Δt)
    nothing
end


"""
    add_baseline_and_extend_tail(wv::RadiationDetectorSignals.RDWaveform{T,U,TV,UV}, n_baseline_samples::Int, total_waveform_length::Int) where {T,U,TV,UV}

Add a zero-valued baseline in front of the waveform `wv` and extends (or cuts off) the waveform at the end with the last value of `wv`.
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
    if length(wv.value) <= total_waveform_length - n_baseline_samples
        new_signal[n_baseline_samples+1:n_baseline_samples+length(wv.value)] = wv.value
        new_signal[n_baseline_samples+length(wv.value)+1:end] .= wv.value[end]
    else
        new_signal[n_baseline_samples+1:end] = wv.value[1:total_waveform_length - n_baseline_samples]
    end
    new_times = if TV <: AbstractRange
        range( zero(first(wv.time)), step = step(wv.time), length = total_waveform_length )
    else
        error("Not yet definted for timestamps of type `$(TV)`")
    end
    return RDWaveform( new_times, new_signal )
end
