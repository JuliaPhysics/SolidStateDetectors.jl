include("table_utils.jl")

"""
    simulate_waveforms( mcevents::TypedTables.Table, sim::Simulation{T})

1. Calculates the drift paths of all energy hits defined in `mcevents`
    based on the drift fields for electrons and holes stored in `s.electron_drift_field` and 
    `s.hole_drift_field`.
2. Determines the signal (waveforms) of all channels 
    (for which a weighting potential is given in the simulation object `s`.)
    for the generated drift paths of the hits.
3. Returns a the input table `mcevents` with an additional column `waveform` 
    in which the generated waveforms are stored. 

Note: The drift paths are just calculated temporarily and not returned. 

# Keywords
- `Δt::RealQuantity = 4u"ns"`: Time difference between two time stamps of the drift and the signals.
- `max_nsteps::Int = 1000`: Maximum number of steps in the drift of each hit. 
- `verbose = false`: Activate or deactivate additional info output. 
"""
function simulate_waveforms( mcevents::TypedTables.Table, sim::Simulation{T};
                             Δt::RealQuantity = 4u"ns",
                             max_nsteps::Int = 1000,
                             verbose = false ) where {T <: SSDFloat}
    n_total_physics_events = length(mcevents)
    Δtime = T(to_internal_units(Δt)) 
    n_contacts = length(sim.detector.contacts)
    S = get_coordinate_system(sim)
    contacts = sim.detector.contacts;
    contact_ids = Int[]
    for contact in contacts if !ismissing(sim.weighting_potentials[contact.id]) push!(contact_ids, contact.id) end end
    wps_interpolated = [ interpolated_scalarfield(sim.weighting_potentials[id]) for id in contact_ids ];
    e_drift_field = get_interpolated_drift_field(sim.electron_drift_field);
    h_drift_field = get_interpolated_drift_field(sim.hole_drift_field);
    
    @info "Detector has $(n_contacts) contact"*(n_contacts != 1 ? "s" : "")
    @info "Table has $(length(mcevents)) physics events ($(sum(map(edeps -> length(edeps), mcevents.edep))) single charge depositions)."

    # First simulate drift paths
    drift_paths = _simulate_charge_drifts(mcevents, sim, Δt, max_nsteps, e_drift_field, h_drift_field, verbose)
    # now iterate over contacts and generate the waveform for each contact
    @info "Generating waveforms..."
    waveforms = map( 
        wp ->  map( 
            x -> _generate_waveform(x.dps, to_internal_units.(x.edeps), Δt, Δtime, wp, S),
            TypedTables.Table(dps = drift_paths, edeps = mcevents.edep)
        ),
        wps_interpolated
    )
    mcevents_chns = map(
        i -> add_column(mcevents, :chnid, fill(contact_ids[i], n_total_physics_events)),
        eachindex(contact_ids)
    )
    mcevents_chns = map(
        i -> add_column(mcevents_chns[i], :waveform, ArrayOfRDWaveforms(waveforms[i])),
        eachindex(waveforms)
    )
    return vcat(mcevents_chns...)  
end


function _simulate_charge_drifts( mcevents::TypedTables.Table, sim::Simulation{T},
                                  Δt::RealQuantity, max_nsteps::Int, 
                                  e_drift_field::Interpolations.Extrapolation, h_drift_field::Interpolations.Extrapolation, 
                                  verbose::Bool ) where {T <: SSDFloat}
    return @showprogress map(mcevents) do phyevt
        _drift_charges(sim.detector, sim.electron_drift_field.grid, sim.point_types, 
                        CartesianPoint{T}.(to_internal_units.(phyevt.pos)),
                        e_drift_field, h_drift_field, 
                        T(Δt.val) * unit(Δt), max_nsteps = max_nsteps, verbose = verbose)
    end
end



function _generate_waveform( drift_paths::Vector{EHDriftPath{T}}, charges::Vector{<:SSDFloat}, Δt::RealQuantity, dt::T,
                             wp::Interpolations.Extrapolation{T, 3}, S::CoordinateSystemType) where {T <: SSDFloat}
    timestamps = _common_timestamps( drift_paths, dt )
    timestamps_with_units = range(zero(Δt), step = Δt, length = length(timestamps))
    signal = zeros(T, length(timestamps))
    add_signal!(signal, timestamps, drift_paths, T.(charges), wp, S)
    RDWaveform( timestamps_with_units, signal )
end


