include("gen_test_charges.jl")
include("geant4_test_charges.jl")
using SolidStateDetectors
import SolidStateDetectors as SSD
using .SSD: drift_charges, _drift_charge!, Event, EHDriftPath, to_internal_units, Event, to_internal_units, interpolated_vectorfield, Electron, internal_length_unit, _is_next_point_in_det, nonzeros
using .SSD: _drift_charge_oldsrp!
using .SSD: ChargeInteractionState, dummy_ci_state, full_ci_state
using .SSD: update_charge_interaction!!, trim_charge_interaction!!, apply_charge_interaction!, resize_charge_interaction_state!!, calc_tmp_d3, _done_field_threshold
using .SSD: ⋮, ⋰, adapted_bcast, adapted_bcast!, parallel_copyto!, parallel_broadcast, parallel_broadcast!
using .SSD: ϵ0, elementary_charge

using LinearAlgebra, SparseArrays
using ArraysOfArrays, StructArrays
using Unitful
using BenchmarkTools
using JLD2

(; edep, pos) = mcevents[1]
pos = CartesianPoint.(to_internal_units.(pos))
charges::Vector{T} = to_internal_units.(edep) ./ to_internal_units(sim.detector.semiconductor.material.E_ionisation)

# For benchmarking _drift_charge!, drift_charges and simulate_waveforms
time_step = 5u"ns"
max_nsteps = 1000
diffusion = false
verbose = true

# self_repulsion = false
#=
self_repulsion = true
=#

srp_trim_threshold = T(20000)
srp_trim_every = 100

# For benchmarking _drift_charge!
electric_field = interpolated_vectorfield(sim.electric_field)
dt::T = T(to_internal_units(time_step))
n_hits::Int = length(pos)
drift_path_e::Array{CartesianPoint{T},2} = Array{CartesianPoint{T},2}(undef, n_hits, max_nsteps)
timestamps_e::Vector{T} = Vector{T}(undef, max_nsteps)
startpos = StructVector(pos)
ϵ_r = T(sim.detector.semiconductor.material.ϵ_r)
current_pos = deepcopy(startpos)

# Benchmark _drift_charge!

_drift_charge!(drift_path_e, timestamps_e, sim.detector, sim.point_types, sim.point_types.grid, current_pos, -charges, dt, electric_field, Electron, diffusion=diffusion, self_repulsion=self_repulsion, verbose=verbose)
_drift_charge_oldsrp!(drift_path_e, timestamps_e, sim.detector, sim.point_types, sim.point_types.grid, current_pos, -charges, dt, electric_field, Electron, diffusion=diffusion, self_repulsion=self_repulsion, verbose=verbose)

@profview @benchmark _drift_charge!(drift_path_e, timestamps_e, sim.detector, sim.point_types, sim.point_types.grid, current_pos, -charges, dt, electric_field, Electron, diffusion=diffusion, self_repulsion=self_repulsion, verbose=verbose, srp_trim_threshold=srp_trim_threshold, srp_trim_every=srp_trim_every)
@profview @benchmark _drift_charge_oldsrp!(drift_path_e, timestamps_e, sim.detector, sim.point_types, sim.point_types.grid, current_pos, -charges, dt, electric_field, Electron, diffusion=diffusion, self_repulsion=self_repulsion, verbose=verbose)

# checking wether all(done) is true and turning self_repulsion off -> this way both drift_charge! are at mean: 89μs (new) and 71μs (old)
# new function only a bit slower because of dummy_ci_state

# Benchmark ChargeInteractionState

field_vectors = fill!(similar(current_pos, CartesianVector{T}), CartesianVector{T}(0, 0, 0))
done = fill!(similar(charges, Bool), false)
field_threshold = srp_trim_threshold


## Dummy adjacency matrix:
ci_state = dummy_ci_state(current_pos, charges)
@profview @benchmark dummy_ci_state($current_pos, $charges)


## Full adjacency matrix:
ci_state = full_ci_state(current_pos, charges, ϵ_r)
ci_state = update_charge_interaction!!(ci_state, current_pos, charges)
apply_charge_interaction!(field_vectors, ci_state)

@benchmark full_ci_state($current_pos, $charges, $ϵ_r)
@benchmark update_charge_interaction!!($ci_state, $current_pos, $charges)
@benchmark apply_charge_interaction!($field_vectors, $ci_state)


## Empty adjacency matrix:
M_adj_empty = sparse(fill!(similar(charges, (length(charges), length(charges))), 0))
ci_state = ChargeInteractionState(current_pos, charges, ϵ_r, M_adj_empty)
ci_state = update_charge_interaction!!(ci_state, current_pos, charges)
ci_state = trim_charge_interaction!!(ci_state, done, srp_trim_threshold) # DEBUG! Produces NaNs!
ci_state = update_charge_interaction!!(ci_state, current_pos, charges)
apply_charge_interaction!(field_vectors, ci_state)

@benchmark update_charge_interaction!!($ci_state, $current_pos, $charges)
@benchmark apply_charge_interaction!($field_vectors, $ci_state)


## Trimmed adjacency matrix:
ci_state = full_ci_state(current_pos, charges, ϵ_r)
ci_state = update_charge_interaction!!(ci_state, current_pos, charges)
ci_state = trim_charge_interaction!!(ci_state, done, srp_trim_threshold)
ci_state = update_charge_interaction!!(ci_state, current_pos, charges)
apply_charge_interaction!(field_vectors, ci_state)

@benchmark update_charge_interaction!!($ci_state, $current_pos, $charges)
@benchmark apply_charge_interaction!($field_vectors, $ci_state)

ci_state = full_ci_state(current_pos, charges, ϵ_r)
ci_state = update_charge_interaction!!(ci_state, current_pos, charges)
@benchmark trim_charge_interaction!!(ci_state, $done, $srp_trim_threshold)

ci_state = full_ci_state(current_pos, charges, ϵ_r)
ci_state = update_charge_interaction!!(ci_state, current_pos, charges)
ci_state = trim_charge_interaction!!(ci_state, done, srp_trim_threshold)
@benchmark update_charge_interaction!!(ci_state, $current_pos, $charges)


# Profiling on all charge interaction steps:
@profview let current_pos = current_pos, charges = charges, ϵ_r = ϵ_r, field_vectors = field_vectors, ci_state = dummy_ci_state(current_pos, charges), done = done, srp_trim_threshold=srp_trim_threshold, t0 = time()
    while time() < t0 + 10
        ci_state = full_ci_state(current_pos, charges, ϵ_r)
        ci_state = update_charge_interaction!!(ci_state, current_pos, charges)
        ci_state = trim_charge_interaction!!(ci_state, done, srp_trim_threshold)
        ci_state = update_charge_interaction!!(ci_state, current_pos, charges)
        apply_charge_interaction!(field_vectors, ci_state)
    end
end


# Benchmark old _add_fieldvector_selfrepulsion!:

n_hits, max_nsteps = size(drift_path_e)
field_vectors = fill!(similar(current_pos, CartesianVector{T}), CartesianVector{T}(0, 0, 0))
SSD._add_fieldvector_selfrepulsion!(field_vectors, current_pos, done, charges, ϵ_r)

@benchmark SSD._add_fieldvector_selfrepulsion!($_field_vectors, $_current_pos, $_done, $_charges, $_ϵ_r)

# old function is fast because all of it is done

# benchmarking drift_charges
n_extra_charges = 1
evt = SolidStateDetectors.Event(pos, edep, n_extra_charges)

@benchmark drift_charges(sim, evt.locations, evt.energies, Δt=time_step, max_nsteps=max_nsteps, diffusion=diffusion, self_repulsion=self_repulsion, verbose=verbose)

#benchmarking simulate_waveforms
# @benchmark simulate_waveforms(mcevents, sim)


# benchmarking _add_fieldvector_selfrepulsion
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
SSD._add_fieldvector_drift!(step_vectors, velocity_vector, current_pos, done, electric_field, det, S)

# Memory allocation for parallelizable _add_fieldvector_selfrepulsion!
velocity_vector = Vector{CartesianVector{T}}(undef, length(step_vectors))

matrix, X, Y, Z = SSD.build_cutoff_matrix(current_pos, electric_field) #, self_repulsion_cutoff

adj_I, adj_J = SSD.findnz(matrix)

# For equivalent of nonzeros(spdiagm(X) * Adj), etc.:
XI = similar(X, length(adj_I)); YI = similar(Y, length(adj_I)); ZI = similar(Z, length(adj_I))
view_XI = view(X, adj_I); view_YI = view(Y, adj_I); view_ZI = view(Z, adj_I)
# For equivalent of nonzeros(Adj * spdiagm(X)), etc.:
XJ = similar(X, length(adj_J)); YJ = similar(Y, length(adj_J)); ZJ = similar(Z, length(adj_J))
view_XJ = view(X, adj_J); view_YJ = view(Y, adj_J); view_ZJ = view(Z, adj_J)

# Δx, Δy, Δz where Adj non-zero:
ΔX_nz = similar(nonzeros(matrix)); ΔY_nz = similar(nonzeros(matrix)); ΔZ_nz = similar(nonzeros(matrix))


Tmp_D3_nz = similar(nonzeros(matrix))

# 1/(4π * ϵ0 * ϵ_r) * Δx/distance^3, same for Δy and Δz:
S_X = similar(matrix); S_Y = similar(matrix); S_Z = similar(matrix)

# E-field components:
Field_X = similar(X); Field_Y = similar(Y); Field_Z = similar(Z)

@benchmark SSD._add_fieldvector_selfrepulsion!(step_vectors, XI, YI, ZI, XJ, YJ, ZJ, view_XI, view_YI, view_ZI, view_XJ, view_YJ, view_ZJ, ΔX_nz, ΔY_nz, ΔZ_nz, Tmp_D3_nz, charges, S_X, S_Y, S_Z, Field_X, Field_Y, Field_Z, ϵ_r)
@profview for i in 1:20
    SSD._add_fieldvector_selfrepulsion!(step_vectors, XI, YI, ZI, XJ, YJ, ZJ, view_XI, view_YI, view_ZI, view_XJ, view_YJ, view_ZJ, ΔX_nz, ΔY_nz, ΔZ_nz, Tmp_D3_nz, charges, S_X, S_Y, S_Z, Field_X, Field_Y, Field_Z, ϵ_r)
end


