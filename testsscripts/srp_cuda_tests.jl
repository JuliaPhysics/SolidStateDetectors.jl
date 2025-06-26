using SolidStateDetectors
import SolidStateDetectors as SSD
using .SSD: to_internal_units, Electron, Hole, internal_length_unit, ϵ0, elementary_charge, getVe, getVh
using .SSD: ChargeInteractionState, dummy_ci_state, full_ci_state
using .SSD: update_charge_interaction!!, trim_charge_interaction!!, apply_charge_interaction!, _get_drift_steps!, _add_fieldvector_selfrepulsion!

using LinearAlgebra, SparseArrays, Random
using ArraysOfArrays, StructArrays, StaticArrays
using Unitful
using BenchmarkTools
using JLD2

using CUDA, Adapt

using ThreadPinning, ParallelProcessingTools
LinearAlgebra.BLAS.set_num_threads(Base.Threads.nthreads())
pinthreads(AutoThreadPinning(random = true, blas = false))

T = Float32
time_step = 4u"ns"
max_nsteps = 100
diffusion = false
self_repulsion = true
verbose = true
self_repulsion_scale = 0.1u"mm"
CC = Electron
n_charges = 10000
edep = 1u"MeV" * randexp(T, n_charges) / n_charges
charge_cloud_scale = 0.1u"mm"
material = :HPGe
cdm = ADLChargeDriftModel()


unitful_startpos = StructVector(CartesianPoint.(
    T.(charge_cloud_scale/2 * randn(n_charges)), T.(charge_cloud_scale/2 * randn(n_charges)), T.(charge_cloud_scale/2 * randn(n_charges))
))

material_properties = SolidStateDetectors.material_properties[material]
cc_sign = CC == Electron ? -1 : 1
E_ionisation = material_properties.E_ionisation
ϵ_r = T(material_properties.ϵ_r)
charges::Vector{T} = cc_sign .* to_internal_units.(edep) ./ to_internal_units(E_ionisation)
startpos = to_internal_units.(StructVector(deepcopy(unitful_startpos)))
current_pos = deepcopy(startpos)

dt::T = T(to_internal_units(time_step))
get_velocity = CC == Electron ? getVe : getVh
steplen_per_Vm = norm(get_velocity(SVector{3,T}(0,0,1), cdm) * dt)
srp_trim_threshold = T(to_internal_units(self_repulsion_scale) / (max_nsteps / 10) / steplen_per_Vm)

field_vectors = fill!(similar(current_pos, CartesianVector{T}), CartesianVector{T}(0, 0, 0))
step_vectors = fill!(similar(current_pos, CartesianVector{T}), CartesianVector{T}(0, 0, 0))
done = fill!(similar(charges, Bool), false)
drift_path_e::Array{CartesianPoint{T},2} = Array{CartesianPoint{T},2}(undef, n_charges, max_nsteps)
timestamps_e::Vector{T} = Vector{T}(undef, max_nsteps)


current_pos = deepcopy(startpos)
ci_state = full_ci_state(current_pos, charges, ϵ_r)
ci_state = update_charge_interaction!!(ci_state, current_pos, charges)
nnz(ci_state.M_adj)

#=
for i in 1:max_nsteps
    t0 = time_ns()
    fill!(field_vectors, zero(eltype(field_vectors)))
    update_charge_interaction!!(ci_state, current_pos, charges)
    apply_charge_interaction!(field_vectors, ci_state)
    ci_state = trim_charge_interaction!!(ci_state, done, srp_trim_threshold)
    _get_drift_steps!(step_vectors, field_vectors, done, dt, cdm, CC)
    current_pos .+= step_vectors
    t1 = time_ns()
end
nnz(ci_state.M_adj)
=#

adapter = CuArray{Float32}
cu_ci_state = adapt(adapter, ci_state)
cu_current_pos = adapt(adapter, current_pos)
cu_charges = adapt(adapter, charges)
cu_field_vectors = adapt(adapter, field_vectors);

update_charge_interaction!!(cu_ci_state, cu_current_pos, cu_charges);
apply_charge_interaction!(cu_field_vectors, cu_ci_state);

# Old charge interaction:
@benchmark begin
    SSD._add_fieldvector_selfrepulsion!($field_vectors, $current_pos, $done, $charges, $ϵ_r)
    sum($field_vectors.x + $field_vectors.y + $field_vectors.z) # Force sync
end

# New charge interaction, CPU:
@benchmark begin
    update_charge_interaction!!($ci_state, $current_pos, $charges)
    apply_charge_interaction!($field_vectors, $ci_state)
    sum($field_vectors.x + $field_vectors.y + $field_vectors.z) # Force sync
end

# New charge interaction, GPU:
@benchmark begin
    update_charge_interaction!!($cu_ci_state, $cu_current_pos, $cu_charges)
    update_charge_interaction!!($cu_ci_state, $cu_current_pos, $cu_charges)
    sum($cu_field_vectors.x + $cu_field_vectors.y + $cu_field_vectors.z) # Force sync
end
