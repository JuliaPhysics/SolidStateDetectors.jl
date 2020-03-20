"""
    mutable struct Event{T <: SSDFloat}

Collection struct for individual events. This (mutable) struct is ment to be used to look at individual events,
not to process a hugh amount of events.
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
    evt.energies = to_internal_units(internal_energy_unit, energies)
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
    if ismissing(T) T = eltype(to_internal_units(internal_energy_unit, evt[:edep][:])) end

    event = Event( 
        CartesianPoint{T}.(to_internal_units.(internal_length_unit, evt[:pos][:])),
        T.(to_internal_units.(internal_energy_unit, evt[:edep][:]))
    )

    return event
end

in(evt::Event, detector::SolidStateDetector) = all( pt -> pt in detector, evt.locations)
in(evt::Event, simulation::Simulation) = all( pt -> pt in simulation.detector, evt.locations)


function drift_charges!(event::Event{T}, sim::Simulation{T}; max_nsteps::Int = 1000, Δt::RealQuantity = 5u"ns", verbose::Bool = true)::Nothing where {T <: SSDFloat}
    event.drift_paths = drift_charges(sim, CartesianPoint.(event.locations), Δt = Δt, max_nsteps = max_nsteps, verbose = verbose)
    nothing
end
function get_signal!(event::Event{T}, sim::Simulation{T}, contact_id::Int; Δt::RealQuantity = 5u"ns")::Nothing where {T <: SSDFloat}
    @assert !ismissing(event.drift_paths) "No drit path for this event. Use `drift_charges(event::Event, sim::Simulation` first."
    @assert !ismissing(sim.weighting_potentials[contact_id]) "No weighting potential for contact $(contact_id). Use `calculate_weighting_potential!(sim, contact_id)` first."
    if ismissing(event.waveforms)
        event.waveforms = Union{Missing, RadiationDetectorSignals.RDWaveform}[missing for i in eachindex(sim.detector.contacts)]
    end
    event.waveforms[contact_id] = get_signal(sim, event.drift_paths, event.energies, contact_id)
    nothing
end
function get_signals!(event::Event{T}, sim::Simulation{T}; Δt::RealQuantity = 5u"ns")::Nothing where {T <: SSDFloat}
    @assert !ismissing(event.drift_paths) "No drit path for this event. Use `drift_charges(event::Event, sim::Simulation` first."
    if ismissing(event.waveforms)
        event.waveforms = Union{Missing, RadiationDetectorSignals.RDWaveform}[missing for i in eachindex(sim.detector.contacts)]
    end
    for contact in sim.detector.contacts
        @assert !ismissing(sim.weighting_potentials[contact.id]) "No weighting potential for contact $(contact.id). Use `calculate_weighting_potential!(sim, contact_id)` first."
        event.waveforms[contact.id] = get_signal(sim, event.drift_paths, event.energies, contact.id, Δt = Δt)
    end
    nothing
end

@recipe function f(wvs::Vector{Union{Missing, RadiationDetectorSignals.RDWaveform}})
    @series begin
        RadiationDetectorSignals.RDWaveform[wv for wv in skipmissing(wvs)]
    end
end

function simulate!(event::Event{T}, sim::Simulation{T}; max_nsteps::Int = 1000, Δt::RealQuantity = 5u"ns", verbose::Bool = true)::Nothing where {T <: SSDFloat}
    drift_charges!(event, sim, max_nsteps = max_nsteps, Δt = Δt, verbose = verbose)
    get_signals!(event, sim, Δt = Δt)
    nothing
end

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