# include("gen_test_charges.jl")
using SolidStateDetectors
using SolidStateDetectors: drift_charges, _drift_charge!, Event, EHDriftPath, to_internal_units, Event, to_internal_units, interpolated_vectorfield, Electron, internal_length_unit,_is_next_point_in_det
using Unitful
using BenchmarkTools
using JLD2


T = Float32
sim = JLD2.load("../test_sim.jld2", "sim")
mcevents = JLD2.load("../test_mcevents.jld2", "mcevents")

(; edep, pos) = mcevents[1]
pos = CartesianPoint.(ustrip.(pos))
charges::Vector{T} = to_internal_units.( edep) ./ to_internal_units(sim.detector.semiconductor.material.E_ionisation)

# benchmarking _drift_charge!, drift_charges and simulate_waveforms
    time_step = 5u"ns"
    max_nsteps = 1000
    diffusion = false
    self_repulsion = true
    verbose = true

        # benchmarking _drift_charge!
            electric_field = interpolated_vectorfield(sim.electric_field)
            dt::T = T(to_internal_units(time_step))
            n_hits::Int = length(pos)
            drift_path_e::Array{CartesianPoint{T},2} = Array{CartesianPoint{T},2}(undef, n_hits, max_nsteps)
            timestamps_e::Vector{T} = Vector{T}(undef, max_nsteps)

            @benchmark _drift_charge!(drift_path_e, timestamps_e, sim.detector, sim.point_types, sim.point_types.grid, pos, -charges, dt, electric_field, Electron, diffusion=diffusion, self_repulsion=self_repulsion, verbose=verbose)
            @profview _drift_charge!(drift_path_e, timestamps_e, sim.detector, sim.point_types, sim.point_types.grid, pos, -charges, dt, electric_field, Electron, diffusion=diffusion, self_repulsion=self_repulsion, verbose=verbose)
            @profview for i in 1:100; _drift_charge!(drift_path_e, timestamps_e, sim.detector, sim.point_types, sim.point_types.grid, pos, -charges, dt, electric_field, Electron, diffusion=diffusion, self_repulsion=self_repulsion, verbose=verbose); end


        # benchmarking drift_charges
            n_extra_charges = 1
            evt = SolidStateDetectors.Event(pos, edep, n_extra_charges)

            @benchmark drift_charges(sim, evt.locations, evt.energies, Δt=time_step, max_nsteps=max_nsteps, diffusion=diffusion, self_repulsion=self_repulsion, verbose=verbose)

        #benchmarking simulate_waveforms
            # @benchmark simulate_waveforms(mcevents, sim)



# benchmarking _add_fieldvector_selfrepulsion!
    det = sim.detector
    point_types = sim.point_types
    charges = -charges
    n_hits::Int, _ = size(drift_path_e)
    ϵ_r::T = T(det.semiconductor.material.ϵ_r)
    S = SSD.Cylindrical
    current_pos::Vector{CartesianPoint{T}} = deepcopy(pos)
    step_vectors::Vector{CartesianVector{T}} = Vector{CartesianVector{T}}(undef, n_hits)
    done::Vector{Bool} = broadcast(pt -> !_is_next_point_in_det(pt, det, point_types), pos)
    SSD._set_to_zero_vector!(step_vectors)
    SSD._add_fieldvector_drift!(step_vectors, current_pos, done, electric_field, det, S)


    @benchmark SSD._add_fieldvector_selfrepulsion!(step_vectors, current_pos, done, charges, ϵ_r)
    @profview for i in 1:100; SSD._add_fieldvector_selfrepulsion!(step_vectors, current_pos, done, charges, ϵ_r); end