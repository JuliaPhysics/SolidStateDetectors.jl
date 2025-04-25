using Plots
using SolidStateDetectors: _check_and_update_position!, get_crossing_pos, project_to_plane, modulate_surface_drift
using Unitful
using Unitful: eV, mm, V, m, C
using PhysicalConstants.CODATA2018: e
using Interpolations
import LegendHDF5IO
using LinearAlgebra
using SparseArrays

include("until_check_and_update_positions.jl")
(; step_vectors, current_pos, done, normal, drift_path, timestamps, istep, det, grid, point_types, startpos, Δ_t, verbose) = g_state
step_vectors



function flatten_vectors(locations)
    return vcat(locations...)
end

function get_cutoff(all_positions; sim=sim)
    arr = []
    N = length(all_positions)
    ϵ₀ = 8.854e-12u"F/m"
    k_e = 1 / (4 * π * ϵ₀)
    el_field_itp = interpolated_vectorfield(sim.electric_field.data, sim.electric_field.grid)
    for i in 1:N, j in 1:N
        if j != i
            q1 = e
            q2 = e
            efield = el_field_itp(all_positions[i]...) * u"V/m"# N/C
            F_el = q1 * efield # N 
            
            x1, y1, z1 = all_positions[i] .* u"mm"
            x2, y2, z2 = all_positions[j] .* u"mm"

            dx2 = (x2 - x1)^2
            dy2 = (y2 - y1)^2
            dz2 = (z2 - z1)^2

            r = sqrt(dx2 + dy2 + dz2)
            r = uconvert(u"m", r)
            F_qq = (k_e * q1 * q2) / r^2

            if norm(F_el) >= F_qq
                nothing
            else
                push!(arr, r)
            end
        end
    end
    return maximum(arr)
end


function build_cutoff_matrix()
    current_pos = flatten_vectors(evt.locations)
    d_cutoff = get_cutoff(current_pos)
    N = length(current_pos)
    row_indices = Int[]
    col_indices = Int[]
    x = Quantity[]
    y = Quantity[]
    z = Quantity[]
    values = Int[]

    for p in current_pos
        push!(x, p[1] * u"mm")
        push!(y, p[2] * u"mm")
        push!(z, p[3] * u"mm")
    end
    for i in 1:N, j in 1:N
        if i != j
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            dz = z[i] - z[j]
            dist = sqrt(dx^2 + dy^2 + dz^2)
            dist = uconvert(u"m", dist)

            if dist < d_cutoff
                push!(row_indices, i)
                push!(col_indices, j)
                push!(values, 1)
            end
        end
    end
    return sparse(row_indices, col_indices, values, N, N)
end

matrix = build_cutoff_matrix()
