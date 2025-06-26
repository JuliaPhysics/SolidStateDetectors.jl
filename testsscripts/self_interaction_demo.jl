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

using Makie, WGLMakie, Observables

T = Float32
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
ci_state_obs = Observable(ci_state)
nnz_history = Observable([nnz(ci_state_obs[].M_adj)])
evaltime_history = Observable([NaN * u"μs"])

fig = let fig = Figure(size = (1600, 800), fontsize = 16)
    pos_obs = map(x -> x.pos, ci_state_obs)
    Pos_x, Pos_y, Pos_z = map(p ->p.x, pos_obs), map(p -> p.y, pos_obs), map(p -> p.z, pos_obs)
    M_adj_obs = map(x -> x.M_adj, ci_state_obs)

    ax1 = Axis3(fig[1:2, 1], title = "Charge positions")
    xlims!(ax1, 3 * minimum(startpos.x), 3 * maximum(startpos.x))
    ylims!(ax1, 3 * minimum(startpos.y), 3 * maximum(startpos.y))
    zlims!(ax1, 3 * minimum(startpos.z), 3 * maximum(startpos.z))
    scatter!(ax1, startpos.x, startpos.y, startpos.z, markersize = 4, alpha = 0.5, color = :blue, label = "start")
    scatter!(ax1, Pos_x, Pos_y, Pos_z, markersize = 4, color = :green, label = "current")

    ax2 = Axis(fig[1:2, 2], aspect = 1, yreversed = true, xaxisposition = :top, title = "Adjacency Matrix")
    heatmap!(ax2, M_adj_obs)

    ax3 = Axis(fig[3, 1], title = "Connections per charge", xlabel = "charges per point", ylabel = "number of connections")
    ylims!(ax3, 0, 1.1 * maximum(vec(sum(M_adj_obs[], dims = 1))))
    scatter!(ax3, -charges, map(M -> vec(sum(M, dims = 1)), M_adj_obs), markersize = 4, color = :red)

    ax3 = Axis(fig[3, 1], title = "Connections per charge", xlabel = "charges per point", ylabel = "number of connections")
    ylims!(ax3, 0, 1.1 * maximum(vec(sum(M_adj_obs[], dims = 1))))
    scatter!(ax3, -charges, map(M -> vec(sum(M, dims = 1)), M_adj_obs), markersize = 4, color = :red)

    ax4 = Axis(fig[3, 2], title = "Non-zero elements in adj. matrix", xlabel = "step", ylabel = "NNZ")
    lines!(ax4, nnz_history)
    xlims!(ax4, 0, max_nsteps)
    ylims!(ax4, 0, 1.1 * nnz(M_adj_obs[]))

    # ax5 = Axis(fig[3, 2], yaxisposition = :right, ylabel = "eval. time", dim2_conversion = Makie.UnitfulConversion(u"ms"; units_in_label=true))
    # xlims!(ax5, 0, max_nsteps)
    # ylims!(ax5, 0u"ms", 1u"ms")
    # hidespines!(ax5)
    # hidexdecorations!(ax5)
    # lines!(ax5, evaltime_history, alpha = 0.3)
    fig
end


function run_sim()
    current_pos .= startpos
    ci_state = full_ci_state(current_pos, charges, ϵ_r)
    ci_state = update_charge_interaction!!(ci_state, current_pos, charges)
    ci_state_obs[] = ci_state
    nnz_history[] = empty!(nnz_history[])
    evaltime_history[] = empty!(evaltime_history[])

    for i in 1:max_nsteps
        t0 = time_ns()
        fill!(field_vectors, zero(eltype(field_vectors)))
        update_charge_interaction!!(ci_state, current_pos, charges)
        apply_charge_interaction!(field_vectors, ci_state)
        ci_state = trim_charge_interaction!!(ci_state, done, srp_trim_threshold)
        _get_drift_steps!(step_vectors, field_vectors, done, dt, cdm, CC)
        current_pos .+= step_vectors
        t1 = time_ns()

        ci_state_obs[] = ci_state
        nnz_history[] = push!(nnz_history[], nnz(ci_state.M_adj))
        evaltime_history[] = push!(evaltime_history[], float(t1 - t0)/100 * u"ns")
        sleep(0.0001)
    end
end


display(fig); run_sim()


# Benchmarks:

#=
@benchmark begin
    update_charge_interaction!!($ci_state, $current_pos, $charges)
    apply_charge_interaction!($field_vectors, $ci_state)
end

@benchmark SSD._add_fieldvector_selfrepulsion!($field_vectors, $current_pos, $done, $charges, $ϵ_r)
=#

#=
(;
    M_adj, adj_row_sums, pos, pos_I, pos_J, pos_vI, pos_vJ, Δpos_IJ,
    charges, charges_J, charges_vJ,
    Tmp_D3_nz, S, ϵ_r, onesT, contr_thresh, contr_thresh_I, contr_thresh_vI
) = ci_state
=#


#=
# CUDA:

using CUDA, Adapt

adapter = CuArray{Float32}
cu_ci_state = adapt(adapter, ci_state)
cu_current_pos = adapt(adapter, current_pos)
cu_charges = adapt(adapter, charges)
cu_field_vectors = adapt(adapter, field_vectors);

update_charge_interaction!!(cu_ci_state, cu_current_pos, cu_charges);

@benchmark update_charge_interaction!!($ci_state, $current_pos, $charges)
@benchmark update_charge_interaction!!($cu_ci_state, $cu_current_pos, $cu_charges)


apply_charge_interaction!(cu_field_vectors, cu_ci_state);

@benchmark apply_charge_interaction!($field_vectors, $ci_state)
@benchmark apply_charge_interaction!($cu_field_vectors, $cu_ci_state)
=#
