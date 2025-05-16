using SolidStateDetectors
using SolidStateDetectors: to_internal_units, get_coordinate_system, interpolated_scalarfield,  interpolated_vectorfield, _convert_internal_energy_to_external_charge, _simulate_charge_drifts, _convertEnergyDepsToChargeDeps, VectorOfArrays, flatview, EHDriftPath, ChargeCarrier, _is_next_point_in_det, RealQuantity, CartesianPoint, SVector
using LegendHDF5IO


abstract type Electron <: ChargeCarrier end 

# region loading simulation data and events
    events_in = lh5open("testsscripts/simulation_output.lh5", "r") do h
        LegendHDF5IO.readdata(h.data_store, "SimulationData")
    end

    include("sim.jl")

    sim.detector = SolidStateDetector(sim.detector, ADLChargeDriftModel(T=T))
    calculate_electric_potential!(sim, refinement_limits = [0.4,0.2,0.1,0.06], verbose = false)
    calculate_electric_field!(sim, n_points_in_φ = 10)
    calculate_weighting_potential!(sim, 1, refinement_limits = [0.4,0.2,0.1,0.06], verbose = false)
# endregion
# region function simulate_waveforms start
    # wf = simulate_waveforms(events_in[1:10], sim, Δt = 1u"ns", max_nsteps = 2000)

    mcevents = events_in[1:10]
    sim = sim 
    Δt = 1u"ns"
    max_nsteps = 2000
    diffusion = false
    self_repulsion  = false
    number_of_carriers = 1
    number_of_shells = 1
    max_interaction_distance = NaN
    verbose = false 

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


# endregion
# region function _simulate_charge_drifts start
        # First simulate drift paths
        #drift_paths_and_edeps = _simulate_charge_drifts(mcevents, sim, Δt, max_nsteps, electric_field, diffusion, self_repulsion, number_of_carriers, number_of_shells, max_interaction_distance, verbose)
    mcevents = mcevents
    sim = sim 
    Δt = Δt
    max_nsteps = max_nsteps
    electric_field = electric_field
    diffusion = diffusion
    self_repulsion = self_repulsion
    number_of_carriers = number_of_carriers
    number_of_shells = number_of_shells
    max_interaction_distance = max_interaction_distance
    verbose = verbose

    # @showprogress map(mcevents) do phyevt
    phyevt = first(mcevents)
    locations, edeps = _convertEnergyDepsToChargeDeps(phyevt.pos, phyevt.edep, sim.detector; number_of_carriers, number_of_shells, max_interaction_distance)
# endregion
# region function _drift_charges start 
    # drift_paths = map( i -> _drift_charges(sim.detector, 
    #                                         sim.electric_field.grid, 
    #                                         sim.point_types, VectorOfArrays(locations[i]), 
    #                                         VectorOfArrays(edeps[i]),
    #                                         electric_field, T(Δt.val) * unit(Δt), 
    #                                         max_nsteps = max_nsteps, 
    #                                         diffusion = diffusion, 
    #                                         self_repulsion = self_repulsion, 
    #                                         verbose = verbose),                 
    #             eachindex(edeps))

    det = sim.detector
    grid = sim.electric_field.grid
    point_types = sim.point_types
    mcevents = mcevents
    electric_field = electric_field
    Δt = T(Δt.val) * unit(Δt)
    max_nsteps = max_nsteps
    diffusion = diffusion
    self_repulsion = self_repulsion
    verbose = verbose
    ############################################### setting the type of starting points and energies#################################################
    all_events = SolidStateDetectors.convert_edepevents(mcevents)

    starting_points = SolidStateDetectors._ustrip_recursive(all_events.Pos)  # change: taking all locations
    energies = SolidStateDetectors._ustrip_recursive(all_events.Q) # change: taking all energies

    drift_paths::Vector{EHDriftPath{T}} = Vector{EHDriftPath{T}}(undef, length(flatview(starting_points)))
    dt::T = T(to_internal_units(Δt))

    drift_path_counter::Int = 0

    start_points::Vector{CartesianPoint{T}} = collect(Iterators.flatten(flatview(starting_points))) # change: taking all starting points

    n_hits::Int = length(start_points)
    charges = flatview(energies) ./ to_internal_units(det.semiconductor.material.E_ionisation)

    drift_path_e::Array{CartesianPoint{T}, 2} = Array{CartesianPoint{T}, 2}(undef, n_hits, max_nsteps)
    drift_path_h::Array{CartesianPoint{T}, 2} = Array{CartesianPoint{T}, 2}(undef, n_hits, max_nsteps)
    timestamps_e::Vector{T} = Vector{T}(undef, max_nsteps)
    timestamps_h::Vector{T} = Vector{T}(undef, max_nsteps)
# end region

# region function _drift_charge! start
    # n_e::Int = _drift_charge!( drift_path_e, timestamps_e, det, point_types, grid, start_points, -charges, dt, electric_field, Electron, diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose )
    # n_h::Int = _drift_charge!( drift_path_h, timestamps_h, det, point_types, grid, start_points,  charges, dt, electric_field, Hole, diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose )

    drift_path::Array{CartesianPoint{T},2} = drift_path_e
    timestamps::Vector{T} = timestamps_e
    det = det
    point_types= point_types
    grid = grid
    startpos::AbstractVector{CartesianPoint{T}} = start_points
    charges = -charges
    Δt = dt
    electric_field = electric_field
    CC = Electron
    diffusion = diffusion
    self_repulsion = self_repulsion
    verbose = verbose
    # returns an Int
                        
    max_nsteps = size(drift_path)
    drift_path[:,1] = startpos
    timestamps[1] = zero(T)
    ϵ_r::T = T(det.semiconductor.material.ϵ_r)
    ###### Diffusion later
    diffusion_length::T = if diffusion
        if CC == Electron && haskey(det.semiconductor.material, :De)
            sqrt(6*_parse_value(T, det.semiconductor.material.De, u"m^2/s") * Δt)
        elseif CC == Hole && haskey(det.semiconductor.material, :Dh)
            sqrt(6*_parse_value(T, det.semiconductor.material.Dh, u"m^2/s") * Δt)
        else 
            @warn "Since v0.9.0, diffusion is modelled via diffusion coefficients `De` (for electrons) and `Dh` (for holes).\n" *
                  "Please update your material properties and pass the diffusion coefficients as `De` and `Dh`.\n" *
                  "You can update it in src/MaterialProperties/MaterialProperties.jl or by overwriting\n" *
                  "`SolidStateDetectors.material_properties` in your julia session and reloading the simulation, e.g.\n
                   SolidStateDetectors.material_properties[:HPGe] = (
                      E_ionisation = 2.95u\"eV\",
                      f_fano = 0.129,
                      ϵ_r = 16.0,
                      ρ = 5.323u\"g*cm^-3\",
                      name = \"High Purity Germanium\",
                      ml = 1.64,
                      mt = 0.0819,
                      De = 200u\"cm^2/s\", # new value 200cm^2/s 
                      Dh = 200u\"cm^2/s\"  # new value 200cm^2/s
                   )\n\n" *
                  "More information can be found at:\n" *
                  "https://juliaphysics.github.io/SolidStateDetectors.jl/stable/man/charge_drift/#Diffusion \n"
            @info "Ignoring diffusion for now"
            diffusion = false
            zero(T)
        end
    else
        zero(T)
    end

    last_real_step_index::Int = 1
    current_pos::Vector{CartesianPoint{T}} = deepcopy(startpos)
    step_vectors::Vector{CartesianVector{T}} = Vector{CartesianVector{T}}(undef, n_hits)
    done::Vector{Bool} = broadcast(pt -> !_is_next_point_in_det(pt, det, point_types), startpos)
    normal::Vector{Bool} = deepcopy(done)
    
    @inbounds for istep in 2:max_nsteps
        last_real_step_index += 1
        _set_to_zero_vector!(step_vectors)
        _add_fieldvector_drift!(step_vectors, current_pos, done, electric_field, det, S)
        self_repulsion && _add_fieldvector_selfrepulsion!(step_vectors, current_pos, done, charges, ϵ_r)
        _get_driftvectors!(step_vectors, done, Δt, det.semiconductor.charge_drift_model, CC)
        diffusion && _add_fieldvector_diffusion!(step_vectors, done, diffusion_length)
        _modulate_driftvectors!(step_vectors, current_pos, det.virtual_drift_volumes)
        _check_and_update_position!(step_vectors, current_pos, done, normal, drift_path, timestamps, istep, det, grid, point_types, startpos, Δt, verbose)
        if all(done) break end
    end

    return last_real_step_index

# end region
# region end of function _drift_charges
    for i in eachindex(start_points)
        drift_paths[drift_path_counter + i] = EHDriftPath{T}( drift_path_e[i,1:n_e], drift_path_h[i,1:n_h], timestamps_e[1:n_e], timestamps_h[1:n_h] )
    end

    drift_path_counter += n_hits

    return drift_paths
# endregion

# region end of function _simulate_charge_drifts
    drift_paths, edeps
    #end

# endregion

# region end of function simulate_waveforms 
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

# endregion
