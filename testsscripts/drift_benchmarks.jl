include("gen_test_charges.jl")

using SolidStateDetectors: drift_charges, drift_charges, Event

(;edep, pos) = mcevents[1]


# ...
_drift_charge!( drift_path_e, timestamps_e, det, point_types, grid, start_points, -charges, dt, electric_field, Electron, diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose ) # ::Int


n_extra_charges = 1
evt = SolidStateDetectors.Event(ustrip.(internal_length_unit, pos), edep, n_extra_charges)
time_step = 5u"ns"
max_nsteps = 1000
drift_charges(sim, evt.locations, evt.energies, Δt = time_step, max_nsteps = max_nsteps, diffusion = false, self_repulsion = true, verbose = true)


#=
simulate_waveforms(mcevents, sim)
=#
