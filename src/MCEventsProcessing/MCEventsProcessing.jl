include("table_utils.jl")

"""
    simulate_waveforms( mcevents::TypedTables.Table, sim::Simulation{T}; kwargs...)
    simulate_waveforms( mcevents::TypedTables.Table, sim::Simulation{T}, output_dir::AbstractString, output_base_name::AbstractString; kwargs...)

Simulates the waveforms for all events defined in `mcevents` for a given [`Simulation`] by

1. calculating the drift paths of all energy hits defined in `mcevents`
    based on the drift fields for electrons and holes stored in `sim.electron_drift_field` and 
    `sim.hole_drift_field`, 
2. determining the signal (waveforms) of all channels, for which a weighting potential 
    is given in the simulation object `sim`, from the generated drift paths of the hits.


## Arguments
* `mcevents::TypedTables.Table`: Table with information about events in the simulated setup.
* `sim::Simulation{T}`: [`Simulation`](@ref) which defines the setup in which the charges in `mcevents` should drift.

If [HDF5.jl](https://github.com/JuliaIO/HDF5.jl) is loaded, this function has additional arguments. 

## Additional Arguments (HDF5)
* `output_dir::AbstractString`: Directory where the HDF5 output file is saved.
* `output_base_name::AbstractString`: Basename of the HDF5 output file, default is `"generated_waveforms"`.

## Keywords
* `max_nsteps::Int = 1000`: Maximum number of steps in the drift of each hit. 
* `Δt::RealQuantity = 4u"ns"`: Time step used for the drift.
* `verbose = false`: Activate or deactivate additional info output.
* `chunk_n_physics_events::Int = 1000` (HDF5 only): Number of events that should be saved in a single HDF5 output file.

## Examples
```julia 
simulate_waveforms(mcevents, sim, Δt = 1u"ns", verbose = false)
# => returns the input table `mcevents` with an additional column `waveform` in which the generated waveforms are stored
```
```julia 
import HDF5
simulate_waveforms(mcevents, sim, "output_dir", "my_basename", Δt = 1u"ns", verbose = false)
# => simulates the charge drift and saves the output to "output_dir/my_basename_evts_xx.h5"
```
!!! note 
    The drift paths are just calculated temporarily and not returned.

!!! note
    Using values with units for `Δt` requires the package [Unitful.jl](https://github.com/PainterQubits/Unitful.jl).
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
    wpots_interpolated = [ interpolated_scalarfield(sim.weighting_potentials[id]) for id in contact_ids ];
    e_drift_field = get_interpolated_drift_field(sim.electron_drift_field);
    h_drift_field = get_interpolated_drift_field(sim.hole_drift_field);
    
    @info "Detector has $(n_contacts) contact"*(n_contacts != 1 ? "s" : "")
    @info "Table has $(length(mcevents)) physics events ($(sum(map(edeps -> length(edeps), mcevents.edep))) single charge depositions)."

    # First simulate drift paths
    drift_paths = _simulate_charge_drifts(mcevents, sim, Δt, max_nsteps, e_drift_field, h_drift_field, verbose)
    # now iterate over contacts and generate the waveform for each contact
    @info "Generating waveforms..."
    waveforms = map( 
        wpot ->  map( 
            x -> _generate_waveform(x.dps, to_internal_units.(x.edeps), Δt, Δtime, wpot, S),
            TypedTables.Table(dps = drift_paths, edeps = mcevents.edep)
        ),
        wpots_interpolated
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
                             wpot::Interpolations.Extrapolation{T, 3}, S::CoordinateSystemType) where {T <: SSDFloat}
    timestamps = _common_timestamps( drift_paths, dt )
    timestamps_with_units = range(zero(Δt), step = Δt, length = length(timestamps))
    signal = zeros(T, length(timestamps))
    add_signal!(signal, timestamps, drift_paths, T.(charges), wpot, S)
    RDWaveform( timestamps_with_units, signal )
end


