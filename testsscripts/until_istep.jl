using Plots
using SolidStateDetectors: drift_charges, _drift_charges, interpolated_vectorfield, to_internal_units, interpolate!, Gridded, Linear, extrapolate, Periodic, EHDriftPath, flatview, _drift_charge!, Electron, Hole, _is_next_point_in_det, _set_to_zero_vector!, _add_fieldvector_drift!, _add_fieldvector_selfrepulsion!, get_coordinate_system, get_velocity_vector, _convert_point, _get_driftvectors!, getVe, SVector, _modulate_driftvectors!, modulate_driftvector, _check_and_update_position!, _parse_value, _add_fieldvector_diffusion!
using Unitful
using Interpolations
import LegendHDF5IO


include("sim.jl")

starting_positions = [CartesianPoint{T}(-0.02, 0.015, 0.04), 
    CartesianPoint{T}(0.015, -0.012, 0.02), 
    CartesianPoint{T}(0.01, -0.025, 0.01)]
energy_depos = T[1460, 609, 1000] * u"keV" # are needed later in the signal generation


evt = Event(starting_positions, energy_depos, 4)

time_step = 5u"ns"

# function: drift_charges!(evt, sim, Δt = time_step)
max_nsteps = 1000
Δt = 5u"ns"
diffusion = true
self_repulsion = false
verbose = true

!in(evt, sim.detector) && move_charges_inside_semiconductor!(evt, sim.detector, verbose=verbose)

## function: drift_charges(sim, evt.locations, evt.energies, Δt = Δt, max_nsteps = max_nsteps, diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose)

starting_positions = evt.locations
energies = evt.energies

### function: _drift_charges(   sim.detector, sim.point_types.grid, sim.point_types, starting_positions, energies, interpolated_vectorfield(sim.electric_field), Δt, max_nsteps = max_nsteps, diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose)
det = sim.detector
grid = sim.point_types.grid
point_types = sim.point_types
starting_points = starting_positions
electric_field = interpolated_vectorfield(sim.electric_field)
drift_paths::Vector{EHDriftPath{T}} = Vector{EHDriftPath{T}}(undef, length(flatview(starting_points)))
dt::T = T(to_internal_units(Δt))

drift_path_counter::Int = 0

# for (i, start_points) in enumerate(starting_points)
# i = 1
start_points = flatview(starting_points)

n_hits::Int = length(start_points)
charges::Vector{T} = flatview(energies) ./ to_internal_units(det.semiconductor.material.E_ionisation)

drift_path_e::Array{CartesianPoint{T},2} = Array{CartesianPoint{T},2}(undef, n_hits, max_nsteps)
drift_path_h::Array{CartesianPoint{T},2} = Array{CartesianPoint{T},2}(undef, n_hits, max_nsteps)
timestamps_e::Vector{T} = Vector{T}(undef, max_nsteps)
timestamps_h::Vector{T} = Vector{T}(undef, max_nsteps)
#### _drift_charge!(drift_path_e, timestamps_e, det, point_types, grid, start_points, -charges, dt, electric_field, Electron, diffusion=diffusion, self_repulsion=self_repulsion, verbose=verbose)
CC = dummy_CC = Electron
drift_path = drift_path_e
timestamps = timestamps_e
startpos = start_points
charges = -charges
Δ_t::T = dt
drift_path[:, 1] = startpos
timestamps[1] = zero(T)
ϵ_r::T = T(det.semiconductor.material.ϵ_r)
diffusion_length::T = if diffusion
    if CC == Electron && haskey(det.semiconductor.material, :De)
        sqrt(6*_parse_value(T, det.semiconductor.material.De, u"m^2/s") * Δ_t)
    elseif CC == Hole && haskey(det.semiconductor.material, :Dh)
        sqrt(6*_parse_value(T, det.semiconductor.material.Dh, u"m^2/s") * Δ_t)
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
S = dummy_S = get_coordinate_system(sim)
last_real_step_index::Int = 1
current_pos::Vector{CartesianPoint{T}} = deepcopy(startpos)
step_vectors::Vector{CartesianVector{T}} = Vector{CartesianVector{T}}(undef, n_hits)
# step_vectors::Matrix{T} = Matrix{T}(undef, n_hits, 3)
done::Vector{Bool} = broadcast(pt -> !_is_next_point_in_det(pt, det, point_types), startpos)
normal::Vector{Bool} = deepcopy(done)

### here 
#@inbounds for istep in 2:max_nsteps
istep = 2 # as dummy
g = Grid
last_real_step_index += 1
# geometry
# _set_to_zero_vector!(step_vectors)
# _add_fieldvector_drift!(step_vectors, current_pos, done, electric_field, det, S)
# self_repulsion && _add_fieldvector_selfrepulsion!(step_vectors, current_pos, done, charges, ϵ_r)
# _get_driftvectors!(step_vectors, done, Δ_t, det.semiconductor.charge_drift_model, CC)
# diffusion && _add_fieldvector_diffusion!(step_vectors, done, diffusion_length)
# _modulate_driftvectors!(step_vectors, current_pos, det.virtual_drift_volumes)

global g_state = (; step_vectors, current_pos, done, normal, electric_field, S, charges, ϵ_r, CC, diffusion_length, drift_path, timestamps, istep, det, grid, point_types, startpos, Δ_t, verbose)

