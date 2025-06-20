include("gen_test_charges.jl")

using SolidStateDetectors: drift_charges, _drift_charge!, Event, EHDriftPath, to_internal_units, Event, to_internal_units, interpolated_vectorfield, Electron

(; edep, pos) = mcevents[1]

# ...
electric_field = interpolated_vectorfield(sim.electric_field)
time_step = 5u"ns"
max_nsteps = 1000

pos = ustrip.(internal_length_unit, pos)
drift_paths::Vector{EHDriftPath{T}} = Vector{EHDriftPath{T}}(undef, length(pos))
dt::T = T(to_internal_units(time_step))

drift_path_counter::Int = 0

n_hits::Int = length(pos)
charges::Vector{T} = to_internal_units.( edep) ./ to_internal_units(sim.detector.semiconductor.material.E_ionisation)

drift_path_e::Array{CartesianPoint{T},2} = Array{CartesianPoint{T},2}(undef, n_hits, max_nsteps)
timestamps_e::Vector{T} = Vector{T}(undef, max_nsteps)
# ...
_drift_charge!(drift_path_e, timestamps_e, sim.detector, sim.point_types, sim.point_types.grid, pos, -charges, dt, electric_field, Electron, diffusion=false, self_repulsion=true, verbose=true)


n_extra_charges = 1
evt = SolidStateDetectors.Event(ustrip.(internal_length_unit, pos), edep, n_extra_charges)

drift_charges(sim, evt.locations, evt.energies, Δt=time_step, max_nsteps=max_nsteps, diffusion=false, self_repulsion=true, verbose=true)


#=
simulate_waveforms(mcevents, sim)
=#
