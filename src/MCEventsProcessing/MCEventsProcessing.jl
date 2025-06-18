include("table_utils.jl")

"""
    simulate_waveforms( mcevents::TypedTables.Table, sim::Simulation{T}; kwargs...)
    simulate_waveforms( mcevents::TypedTables.Table, sim::Simulation{T}, output_dir::AbstractString, output_base_name::AbstractString; kwargs...)

Simulates the waveforms for all events defined in `mcevents` for a given [`Simulation`](@ref) by

1. calculating the drift paths of all energy hits defined in `mcevents`,
2. determining the signal (waveforms) for each [`Contact`](@ref), 
    for which a [`WeightingPotential`](@ref) is specified in `sim.weighting_potentials`.

## Arguments
* `mcevents::TypedTables.Table`: Table with information about events in the simulated setup.
* `sim::Simulation{T}`: [`Simulation`](@ref) which defines the setup in which the charges in `mcevents` should drift.

If [LegendHDF5IO.jl](https://github.com/legend-exp/LegendHDF5IO.jl) is loaded, this function has additional arguments. 

## Additional Arguments (`LegendHDF5IO`)
* `output_dir::AbstractString`: Directory where the HDF5 output file is saved.
* `output_base_name::AbstractString`: Basename of the HDF5 output file, default is `"generated_waveforms"`.

## Keywords
* `max_nsteps::Int = 1000`: Maximum number of steps in the drift of each hit. 
* `Δt::RealQuantity = 4u"ns"`: Time step used for the drift.
* `diffusion::Bool = false`: Activate or deactive diffusion of charge carriers via random walk.
* `self_repulsion::Bool = false`: Activate or deactive self-repulsion of charge carriers of the same type.
* `number_of_carriers::Int = 1`: Number of charge carriers to be used in the N-Body simulation of an energy deposition. 
* `number_of_shells::Int = 1`: Number of shells around the `center` point of the energy deposition.
* `max_interaction_distance = NaN`: Maximum distance for which charge clouds will interact with each other (if `NaN`, then all charge clouds drift independently).
* `verbose = false`: Activate or deactivate additional info output.
* `chunk_n_physics_events::Int = 1000` (`LegendHDF5IO` only): Number of events that should be saved in a single HDF5 output file.

## Examples
```julia 
simulate_waveforms(mcevents, sim, Δt = 1u"ns", verbose = false)
# => returns the input table `mcevents` with an additional column `waveform` in which the generated waveforms are stored
```
```julia 
using LegendHDF5IO
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
                             diffusion::Bool = false,
                             self_repulsion::Bool = false,
                             number_of_carriers::Int = 1,
                             number_of_shells::Int = 1,
                             max_interaction_distance::Union{<:Real, <:LengthQuantity} = NaN,
                             verbose::Bool = false ) where {T <: SSDFloat}
    n_total_physics_events = length(mcevents)
    Δtime = T(to_internal_units(Δt)) 
    n_contacts = length(sim.detector.contacts)
    S = get_coordinate_system(sim)
    contacts = sim.detector.contacts;
    contact_ids = Int[]
    for contact in contacts if !ismissing(sim.weighting_potentials[contact.id]) push!(contact_ids, contact.id) end end
    wpots_interpolated = [ interpolated_scalarfield(sim.weighting_potentials[id]) for id in contact_ids ];
    electric_field = interpolated_vectorfield(sim.electric_field)
    ctm = sim.detector.semiconductor.charge_trapping_model

    unitless_energy_to_charge = _convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)
    @info "Detector has $(n_contacts) contact"*(n_contacts != 1 ? "s" : "")
    @info "Table has $(length(mcevents)) physics events ($(sum(map(edeps -> length(edeps), mcevents.edep))) single charge depositions)."

    # First simulate drift paths
    drift_paths_and_edeps = _simulate_charge_drifts(mcevents, sim, Δt, max_nsteps, electric_field, 
        diffusion, self_repulsion, number_of_carriers, number_of_shells, max_interaction_distance, verbose)
    drift_paths = map(x -> vcat([vcat(ed...) for ed in x[1]]...), drift_paths_and_edeps)
    edeps = map(x -> vcat([vcat(ed...) for ed in x[2]]...), drift_paths_and_edeps)
    # now iterate over contacts and generate the waveform for each contact
    @info "Generating waveforms..."
    waveforms = map( 
        wpot ->  map( 
            x -> _generate_waveform(x.dps, to_internal_units.(x.edeps), Δt, Δtime, wpot, S, unitless_energy_to_charge, ctm),
            TypedTables.Table(dps = drift_paths, edeps = edeps)
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


# Support detectors hits with positions given as vectors:
function _get_unitless_positions(point_vectors::AbstractVector{<:StaticVector{3}})
    Ref(cartesian_zero) .+ to_internal_units.(point_vectors)
end

# Support detectors hits with positions given as points:
_get_unitless_positions(points::AbstractVector{<:CartesianPoint}) = to_internal_units.(points)

function _convertEnergyDepsToChargeDeps(
        positions::AbstractVector, edep::AbstractVector{<:Quantity}, det::SolidStateDetector{T}; 
        particle_type::Type{PT} = Gamma, 
        radius::Vector{<:Union{<:Real, <:LengthQuantity}} = radius_guess.(T.(to_internal_units.(edep)), particle_type),
        max_interaction_distance::Union{<:Real, <:LengthQuantity} = NaN,
        kwargs...
    ) where {T <: SSDFloat, PT <: ParticleType}

    @assert eachindex(positions) == eachindex(edep) == eachindex(radius)

    d::T = T(to_internal_units(max_interaction_distance))
    @assert isnan(d) || d >= 0 "Max. interaction distance must be positive or NaN (no grouping), but $(max_interaction_distance) was given."

    unitless_pos =_get_unitless_positions(positions)
    unitless_edep = T.(to_internal_units.(edep))
    unitless_radius = T.(to_internal_units.(radius))
    loc, ene, rad = group_points_by_distance(unitless_pos, unitless_edep, unitless_radius, d)
    _convertEnergyDepsToChargeDeps(loc, ene, det; radius = rad, kwargs...)
end

function _convertEnergyDepsToChargeDeps(
        pos::AbstractVector{<:AbstractVector{<:CartesianPoint}}, edep::AbstractVector{<:AbstractVector}, det::SolidStateDetector{T}; 
        particle_type::Type{PT} = Gamma,
        radius::AbstractVector{<:AbstractVector{<:Union{<:Real, <:LengthQuantity}}} = map(e -> radius_guess.(to_internal_units.(e), particle_type), edep),
        number_of_carriers::Int = 1, number_of_shells::Int = 1, 
        max_interaction_distance::Union{<:Real, <:LengthQuantity} = NaN # is ignored here
    ) where {T <: SSDFloat, PT <: ParticleType}

    charge_clouds = broadcast(
        iEdep_indep -> broadcast(
            i_together -> begin 
                nbcc = NBodyChargeCloud(
                    CartesianPoint{T}(to_internal_units(pos[iEdep_indep][i_together])), 
                    T.(to_internal_units(edep[iEdep_indep][i_together])),
                    number_of_carriers,
                    number_of_shells = number_of_shells,
                    radius = T(to_internal_units(radius[iEdep_indep][i_together]))
                )
                move_charges_inside_semiconductor!([nbcc.locations], [nbcc.energies], det)
                nbcc
            end,
            eachindex(edep[iEdep_indep])
        ), 
        eachindex(edep)
    )
    locations = map.(cc -> cc.locations, charge_clouds)
    edeps = map.(cc -> cc.energies, charge_clouds)
    locations, edeps
end

function _simulate_charge_drifts( mcevents::TypedTables.Table, sim::Simulation{T},
                                  Δt::RealQuantity, max_nsteps::Int, 
                                  electric_field::Interpolations.Extrapolation,
                                  diffusion::Bool,
                                  self_repulsion::Bool,
                                  number_of_carriers::Int, 
                                  number_of_shells::Int,
                                  max_interaction_distance::Union{<:Real, <:LengthQuantity},
                                  verbose::Bool ) where {T <: SSDFloat}
    @showprogress map(mcevents) do phyevt
        locations, edeps = _convertEnergyDepsToChargeDeps(phyevt.pos, phyevt.edep, sim.detector; number_of_carriers, number_of_shells, max_interaction_distance)
        drift_paths = map( i -> _drift_charges(sim.detector, sim.electric_field.grid, sim.point_types, 
                VectorOfArrays(locations[i]), VectorOfArrays(edeps[i]),
                electric_field, T(Δt.val) * unit(Δt), max_nsteps = max_nsteps, 
                diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose
            ),
            eachindex(edeps)            
        )
        drift_paths, edeps
    end
end



function _generate_waveform( drift_paths::Vector{<:EHDriftPath{T}}, charges::Vector{<:SSDFloat}, Δt::RealQuantity, dt::T,
                             wpot::Interpolations.Extrapolation{T, 3}, S::CoordinateSystemType, unitless_energy_to_charge, 
                             ctm::AbstractChargeTrappingModel{T} = NoChargeTrappingModel{T}()) where {T <: SSDFloat}
    timestamps = _common_timestamps( drift_paths, dt )
    timestamps_with_units = range(zero(Δt), step = Δt, length = length(timestamps))
    signal = zeros(T, length(timestamps))
    add_signal!(signal, timestamps, drift_paths, T.(charges), wpot, S, ctm)
    RDWaveform( timestamps_with_units, signal * unitless_energy_to_charge)
end


"""
    run_geant4_simulation(app::Geant4.G4JLApplication, number_of_events::Int; energy_threshold::Unitful.Energy)

Simulates the given `Geant4.G4JLApplication` until `number_of_events` events with non-zero energy depositions have been generated.

## Arguments
* `app::Geant4.G4JLApplication`: Contains information about the detector setup and the particle source.
* `number_of_events::Int`: The amount of events inside of the detector that should be recorded.

## Keywords
* `energy_threshold::Unitful.Energy`: Defines a lower threshold on the summed energy of an event to be recorded.
   Default is `eps(Float64)*u"keV"`, i.e. all events with non-zero deposits are recorded.
   If set to `0u"keV"`, all primary events will be recorded.

## Examples
```julia 
run_geant4_simulation(app, 1000)
# => returns all events in a `TypedTables.Table` with the fields `evtno`, `detno`, `thit, `edep` and `pos`
```

!!! note 
    Since the function is running until enough events have been collected, setups in which the source
        does not irradiate the detector will run indefinitely. 
"""
function run_geant4_simulation end