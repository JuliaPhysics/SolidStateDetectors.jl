using SolidStateDetectors
import SolidStateDetectors as SSD
using .SSD: drift_charges, _drift_charge!, Event, EHDriftPath, to_internal_units, Event, to_internal_units, interpolated_vectorfield, Electron, internal_length_unit, _is_next_point_in_det, nonzeros
using .SSD: _drift_charge_oldsrp!
using .SSD: ChargeInteractionState, dummy_ci_state, full_ci_state, getVe, getVh
using .SSD: update_charge_interaction!!, trim_charge_interaction!!, apply_charge_interaction!, _get_drift_steps!, _add_fieldvector_selfrepulsion!
using .SSD: resize_charge_interaction_state!!, calc_tmp_d3
using .SSD: ⋮, ⋰, adapted_bcast, adapted_bcast!, parallel_copyto!, parallel_broadcast, parallel_broadcast!
using .SSD: ϵ0, elementary_charge

using LinearAlgebra, SparseArrays, Random
using ArraysOfArrays, StructArrays
using Unitful
using BenchmarkTools
using JLD2

using Plots

include("snapshotting.jl")

T = Float32


function build_detsim(::Type{T}, example_geometry::Symbol) where T<:Real
    detector_config_filename = SSD_examples[example_geometry]
    sim = Simulation{T}(detector_config_filename)
    calculate_electric_potential!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1, 0.05, 0.01])
    calculate_electric_field!(sim, n_points_in_φ = 72)
    for contact in sim.detector.contacts
        calculate_weighting_potential!(sim, contact.id, refinement_limits = [0.2, 0.1, 0.05, 0.01], n_points_in_φ = 2, verbose = false)
    end
    charge_drift_model = ADLChargeDriftModel()
    sim.detector = SolidStateDetector(sim.detector, charge_drift_model);
    return sim
end

@snapshot sim = build_detsim(Float32, :InvertedCoax)


time_step = 4u"ns"
max_nsteps = 1000
diffusion = false
self_repulsion = true
verbose = true
self_repulsion_scale = 0.1u"mm"
CC = Electron
n_charges = 1000
edep = 1u"MeV" * randexp(T, n_charges) / n_charges
charge_cloud_scale = 0.1u"mm"


unitful_startpos = StructVector(CartesianPoint.(
    T.(charge_cloud_scale * randn(n_charges)), T.(charge_cloud_scale * randn(n_charges)), T.(charge_cloud_scale * randn(n_charges))
))

cc_sign = CC == Electron ? -1 : 1
charges::Vector{T} = cc_sign .* to_internal_units.(edep) ./ to_internal_units(sim.detector.semiconductor.material.E_ionisation)
startpos = to_internal_units.(StructVector(deepcopy(unitful_startpos)))
current_pos = deepcopy(startpos)

electric_field = interpolated_vectorfield(sim.electric_field)
dt::T = T(to_internal_units(time_step))
ϵ_r = T(sim.detector.semiconductor.material.ϵ_r)
cdm = sim.detector.semiconductor.charge_drift_model
steplen_per_Vm = norm(getVe(SVector{3,T}(0,0,1), cdm) * dt)
srp_trim_threshold = to_internal_units(self_repulsion_scale) / (max_nsteps / 10) / steplen_per_Vm

field_vectors = fill!(similar(current_pos, CartesianVector{T}), CartesianVector{T}(0, 0, 0))
step_vectors = fill!(similar(current_pos, CartesianVector{T}), CartesianVector{T}(0, 0, 0))
done = fill!(similar(charges, Bool), false)
drift_path_e::Array{CartesianPoint{T},2} = Array{CartesianPoint{T},2}(undef, n_charges, max_nsteps)
timestamps_e::Vector{T} = Vector{T}(undef, max_nsteps)




# _drift_charge!(drift_path_e, timestamps_e, sim.detector, sim.point_types, sim.point_types.grid, current_pos, -charges, dt, electric_field, CC, diffusion=diffusion, self_repulsion=self_repulsion, verbose=verbose)

ci_state = full_ci_state(current_pos, charges, ϵ_r)
current_pos = deepcopy(startpos)
nonzeros(ci_state.M_adj)


ci_state = update_charge_interaction!!(ci_state, current_pos, charges)
fill!(field_vectors, zero(eltype(field_vectors)))
apply_charge_interaction!(field_vectors, ci_state)
ci_state = trim_charge_interaction!!(ci_state, done, srp_trim_threshold)
_get_drift_steps!(step_vectors, field_vectors, done, dt, cdm, CC)
current_pos .+= step_vectors
length(nonzeros(ci_state.M_adj))

#scatter((x -> x.x).(startpos), (x -> x.y).(startpos)); scatter!((x -> x.x).(current_pos), (x -> x.y).(current_pos))
#stephist(norm.(startpos .- Ref(cartesian_zero)) .* 1000, bins = 0:0.01:1, label = "start", xlabel = "distance from origin / mm")
#stephist!(norm.(current_pos .- Ref(cartesian_zero)) .* 1000, bins = 0:0.01:1, label = "current")

plot(startpos); plot!(current_pos)


(;
    M_adj, adj_row_sums, pos, pos_I, pos_J, pos_vI, pos_vJ, Δpos_IJ,
    charges, charges_J, charges_vJ,
    Tmp_D3_nz, S, ϵ_r, onesT, contr_thresh, contr_thresh_I, contr_thresh_vI
) = ci_state
