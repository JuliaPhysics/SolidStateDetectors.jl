"""
    mutable struct Event{T <: SSDFloat}

Collection struct for individual events. This (mutable) struct is ment to be used to look at individual events,
not to process a hugh amount of events.
"""
mutable struct Event{T <: SSDFloat} 
    locations::Union{Vector{<:AbstractCoordinatePoint{T}}, Missing}
    energies::Union{Vector{T}, Missing}
    drift_paths::Union{Vector{DriftPath{T}}, Missing}
    signals::Union{Missing, Vector{Any}}
    Event{T}() where {T <: SSDFloat} = new{T}(missing, missing, missing, missing)
end

function Event(locations::Vector{<:AbstractCoordinatePoint{T}}, energies::Vector{<:RealQuantity{T}} = ones(T, length(locations)))::Event{T} where {T <: SSDFloat}
    evt = Event{T}()
    evt.locations = locations
    evt.energies = to_internal_units(internal_energy_unit, energies)
    evt.signals = missing
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

    evt = Event( 
        CartesianPoint{T}.(to_internal_units.(internal_length_unit, evt[:pos][:])),
        T.(to_internal_units.(internal_energy_unit, evt[:edep][:]))
    )

    return evt
end

function drift_charges!(event::Event{T}, setup::SSDSetup{T}; n_steps::Int = 1000, Δt::RealQuantity = 5u"ns", verbose::Bool = true)::Nothing where {T <: SSDFloat}
    event.drift_paths = drift_charges(setup, CartesianPoint.(event.locations), Δt = Δt, n_steps = n_steps, verbose = verbose)
    nothing
end
function get_signal!(event::Event{T}, setup::SSDSetup{T}, contact_id::Int)::Nothing where {T <: SSDFloat}
    @assert !ismissing(event.drift_paths) "No drit path for this event. Use `drift_charges(event::Event, setup::SSDSetup` first."
    @assert !ismissing(setup.weighting_potentials[contact_id]) "No weighting potential for contact $(contact_id). Use `calculate_weighting_potential!(setup, contact_id)` first."
    if ismissing(event.signals)
        event.signals = Any[missing for i in eachindex(setup.detector.contacts)]
    end
    event.signals[contact_id] = get_signal(setup, event.drift_paths, event.energies, contact_id)
    nothing
end
function get_signals!(event::Event{T}, setup::SSDSetup{T})::Nothing where {T <: SSDFloat}
    @assert !ismissing(event.drift_paths) "No drit path for this event. Use `drift_charges(event::Event, setup::SSDSetup` first."
    if ismissing(event.signals)
        event.signals = Any[missing for i in eachindex(setup.detector.contacts)]
    end
    for contact in setup.detector.contacts
        @assert !ismissing(setup.weighting_potentials[contact.id]) "No weighting potential for contact $(contact.id). Use `calculate_weighting_potential!(setup, contact_id)` first."
        event.signals[contact.id] = get_signal(setup, event.drift_paths, event.energies, contact.id)
    end
    nothing
end

function simulate!(event::Event{T}, setup::SSDSetup{T}; n_steps::Int = 1000, Δt::RealQuantity = 5u"ns", verbose::Bool = true)::Nothing where {T <: SSDFloat}
    drift_charges!(event, setup, n_steps = n_steps, Δt = Δt, verbose = verbose)
    get_signals!(event, setup)
    nothing
end