using Plots
using SolidStateDetectors: elementary_charge, ϵ0, distance_squared, normalize, StaticVector, _add_fieldvector_selfrepulsion!, SSDFloat, CartesianVector, CartesianPoint, _set_to_zero_vector!, getVe, modulate_driftvector, _get_driftvectors!, get_velocity_vector, _convert_point, _check_and_update_position!, _parse_value, _add_fieldvector_diffusion!, interpolated_vectorfield
using Unitful
using Unitful: eV, mm, V, m, C
using PhysicalConstants.CODATA2018: e
using Interpolations
using LinearAlgebra
using SparseArrays
using StaticArrays, StructArrays, ArraysOfArrays
using Test 
import LegendHDF5IO
import SolidStateDetectors: _add_fieldvector_drift!, _add_fieldvector_selfrepulsion! 
using BenchmarkTools


const CD_ELECTRODE = 0x00
const CD_OUTSIDE = 0x01
const CD_BULK = 0x02
const CD_FLOATING_BOUNDARY = 0x04
const k = Float32(8.9875517923e9)  # N·m²/C² for example

include("until_istep.jl")
(; step_vectors, current_pos, done, normal, electric_field, S, charges, ϵ_r, CC, diffusion_length, drift_path, timestamps, istep, det, grid, point_types, startpos, Δ_t, verbose) = g_state



# Simplified version of SolidStateDetectors.CartesianVector
# struct CartesianVector{T} <: StaticArrays.FieldVector{3, T}
#     x::T
#     y::T
#     z::T
# end

T = Float32

n_events = 100
n_charges_avg = 1000

EvtNo = Vector{Int}()
Q = Vector{T}()
Pos = StructVector{CartesianVector{T}}(undef, 0)

for evtno in 1:n_events
    n_charges = rand(round(Int, 0.5 * n_charges_avg):round(Int, 1.5 * n_charges_avg))
    pos_center = CartesianVector(10 * rand(T), 10 * rand(T), 10 * rand(T))

    Q_evt = 1000 * rand(T, n_charges)
    Pos_evt = [CartesianVector(pos_center.x + 0.2*randn(T), pos_center.y + 0.2*randn(T), pos_center.z + 0.2*randn(T)) for _ in 1:n_charges]
    append!(EvtNo, fill(evtno, n_charges))
    append!(Q, Q_evt)
    append!(Pos, Pos_evt)
end


### Oli version 
function build_adj_matrix(T, EvtNo, Pos, interact_dist)
    n = length(EvtNo)
    nz_I = Int[]
    nz_J = Int[]

    nested_idxs = consgroupedview(EvtNo, eachindex(EvtNo))
    for evt_idxs in nested_idxs
        for i in evt_idxs
            for j in evt_idxs
                if i != j
                    dist = norm(Pos[i] - Pos[j])
                    if dist < interact_dist
                        push!(nz_I, i)
                        push!(nz_J, j)
                    end
                end
            end
        end    
    end
    Adj = sparse(nz_I, nz_J, fill(one(T), length(nz_I)), n, n)
    return Adj
end
# my version
function build_cutoff_matrix(current_pos, d_cutoff)
    #d_cutoff = get_cutoff(current_pos)
    N = length(current_pos)
    row_indices = Int[]
    col_indices = Int[]
    x = T[]
    y = T[]
    z = T[]
    values = Bool[]

    for p in current_pos
        push!(x, p[1])#* u"mm"
        push!(y, p[2])#* u"mm"
        push!(z, p[3])#* u"mm"
    end
    for i in 1:N, j in 1:N
        if i != j
            Δx = x[i] - x[j]
            Δy = y[i] - y[j]
            Δz = z[i] - z[j]
            dist = sqrt(Δx^2 + Δy^2 + Δz^2)/1000

            if dist < d_cutoff
                push!(row_indices, i)
                push!(col_indices, j)
                push!(values, true)
            end
        end
    end
    return sparse(row_indices, col_indices, fill(one(T), length(row_indices)), N, N), x, y, z
end

interact_dist = 0.3
# Adj = build_adj_matrix(T, EvtNo, Pos, interact_dist)

build_adj_matrix(T, EvtNo, Pos, interact_dist)
@time build_cutoff_matrix(Pos, interact_dist)

