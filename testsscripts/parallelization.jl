
# using Plots
using SolidStateDetectors: elementary_charge, ϵ0, distance_squared, normalize, StaticVector, _add_fieldvector_selfrepulsion!, SSDFloat, CartesianVector, CartesianPoint, _set_to_zero_vector!, getVe, modulate_driftvector, _get_driftvectors!, get_velocity_vector, _convert_point, _check_and_update_position!, _parse_value, _add_fieldvector_diffusion!, interpolated_vectorfield
using Unitful
using Unitful: eV, mm, V, m, C
using PhysicalConstants.CODATA2018: e
using Interpolations
using LinearAlgebra
using SparseArrays 
import LegendHDF5IO
import SolidStateDetectors: _add_fieldvector_drift!, _add_fieldvector_selfrepulsion! 

include("until_istep.jl")
(; step_vectors, current_pos, done, normal, electric_field, S, charges, ϵ_r, CC, diffusion_length, drift_path, timestamps, istep, det, grid, point_types, startpos, Δ_t, verbose) = g_state

# _add_fieldvector_drift!(step_vectors, current_pos, done, electric_field, det, S)
# self_repulsion && _add_fieldvector_selfrepulsion!(step_vectors, current_pos, done, charges, ϵ_r)
# _get_driftvectors!(step_vectors, done, Δ_t, det.semiconductor.charge_drift_model, CC)
# diffusion && _add_fieldvector_diffusion!(step_vectors, done, diffusion_length)
# _modulate_driftvectors!(step_vectors, current_pos, det.virtual_drift_volumes)


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
            efield = el_field_itp(all_positions[i]...) * u"V/m"
            F_el = q1 * efield 
            
            x1, y1, z1 = all_positions[i] .* u"mm"
            x2, y2, z2 = all_positions[j] .* u"mm"

            Δx2 = (x2 - x1)^2
            Δy2 = (y2 - y1)^2
            Δz2 = (z2 - z1)^2

            r = sqrt(Δx2 + Δy2 + Δz2)
            r = uconvert(u"m", r)
            F_qq = (k_e * q1 * q2) / r^2

            if norm(F_el) >= F_qq
                nothing
            else
                push!(arr, r)
            end
        end
    end
    return ustrip(maximum(arr))
end

function build_cutoff_matrix(current_pos)
    d_cutoff = get_cutoff(current_pos)
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

function update_distances(XI, YI, ZI, XJ, YJ, ZJ, ΔX_nz, ΔY_nz, ΔZ_nz)
    ΔX_nz .= XI .- XJ
    ΔY_nz .= YI .- YJ
    ΔZ_nz .= ZI .- ZJ
end

# parallelization finished
function _add_fieldvector_drift!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, done::Vector{Bool}, electric_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, det::SolidStateDetector{T}, ::Type{S})::Nothing where {T, S}
    velocity_vector = Vector{CartesianVector{Float32}}(undef, length(step_vectors))# defining before istep
    _set_to_zero_vector!(velocity_vector)
    velocity_vector .= get_velocity_vector.(Ref(electric_field), _convert_point.(current_pos, S))
    step_vectors .= step_vectors .+ velocity_vector .* (.! done)
    nothing
end
# parallelization finished
function _add_fieldvector_selfrepulsion!(step_vectors::Vector{CartesianVector{T}}, XI::Vector{T}, YI::Vector{T}, ZI::Vector{T}, XJ::Vector{T}, YJ::Vector{T}, ZJ::Vector{T}, view_XI, view_YI, view_ZI, view_XJ, view_YJ, view_ZJ, ΔX_nz::Vector{T}, ΔY_nz::Vector{T}, ΔZ_nz::Vector{T}, charges::Vector{T}, S_X, S_Y, S_Z, Field_X::Vector{T}, Field_Y::Vector{T}, Field_Z::Vector{T}, ϵ_r::T)::Nothing where {T <: SSDFloat}

    # Update equivalents of nonzeros(spdiagm(X) * Adj), etc.:
    XI .= view_XI
    YI .= view_YI
    ZI .= view_ZI
    # Update equivalents of nonzeros(Adj * spdiagm(X)), etc.:
    XJ .= view_XJ
    YJ .= view_YJ
    ZJ .= view_ZJ

    update_distances(XI, YI, ZI, XJ, YJ, ZJ, ΔX_nz, ΔY_nz, ΔZ_nz)
    
    inv_4π_ϵ0_ϵ_r = inv(4π * ϵ0 * ϵ_r)

    function calc_tmp_d3(Δx::T, Δy::T, Δz::T) where {T<:Real}
        #inv_distance_3 = inv(max(sqrt(Δx * Δx + Δy * Δy + Δz * Δz)^3, T(1e-10)))
        inv_distance_3 = inv(max((Δx^2 + Δy^2 + Δz^2)^(3/2), T(1e-10)))
        elementary_charge * inv_4π_ϵ0_ϵ_r * inv_distance_3
    end

    # Update Tmp_D3_nz:
    Tmp_D3_nz .= calc_tmp_d3.(ΔX_nz, ΔY_nz, ΔZ_nz)
    # Update S_X, S_Y, S_Z:
    nonzeros(S_X) .= ΔX_nz .* Tmp_D3_nz
    nonzeros(S_Y) .= ΔY_nz .* Tmp_D3_nz
    nonzeros(S_Z) .= ΔZ_nz .* Tmp_D3_nz

    # Update E-field components:
    mul!(Field_X, S_X, charges)
    mul!(Field_Y, S_Y, charges)
    mul!(Field_Z, S_Z, charges)

    step_vectors .+= CartesianVector.(Field_X, Field_Y, Field_Z)
    nothing
end

###################################### TODOS
# function _get_driftvectors!(step_vectors::Vector{CartesianVector{T}}, done::Vector{Bool}, Δt::T, cdm::AbstractChargeDriftModel{T}, ::Type{Electron})::Nothing where {T <: SSDFloat}
# cdm = det.semiconductor.charge_drift_model
#     for n in eachindex(step_vectors)
#         if !done[n]
#             step_vectors[n] = getVe(SVector{3,T}(step_vectors[n]), cdm) * Δt
#         end
#     end
#     nothing
# end

# function _add_fieldvector_diffusion!(step_vectors::Vector{CartesianVector{T}}, done::Vector{Bool}, length::T = T(0.5e3))::Nothing where {T <: SSDFloat}
# sinθ::T, cosθ::T = sincos(acos(T(2*rand() - 1)))
# sinφ::T, cosφ::T = sincos(T(rand()*2π))
# dif_vector = CartesianVector{T}( diffusion_length * cosφ * sinθ, diffusion_length * sinφ * sinθ, diffusion_length * cosθ )
# step_vectors .= step_vectors .+ dif_vector' .* (.! done)
# nothing 
# end

# function _modulate_driftvectors!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, vdv::Vector{V})::Nothing where {T <: SSDFloat, V <: AbstractVirtualVolume{T}}
#     for n in eachindex(step_vectors)
#         step_vectors[n] = modulate_driftvector(step_vectors[n], current_pos[n], vdv)
#     end
#     nothing
# end 
###################################### TODOS

velocity_vector = Vector{CartesianVector{T}}(undef, length(step_vectors))

matrix, X, Y, Z = build_cutoff_matrix(current_pos)

nz_I, nz_J = findnz(matrix)

# For equivalent of nonzeros(spdiagm(X) * Adj), etc.:
XI = similar(X, length(nz_I)); YI = similar(Y, length(nz_I)); ZI = similar(Z, length(nz_I))
view_XI = view(X, nz_I); view_YI = view(Y, nz_I); view_ZI = view(Z, nz_I)
# For equivalent of nonzeros(Adj * spdiagm(X)), etc.:
XJ = similar(X, length(nz_J)); YJ = similar(Y, length(nz_J)); ZJ = similar(Z, length(nz_J))
view_XJ = view(X, nz_J); view_YJ = view(Y, nz_J); view_ZJ = view(Z, nz_J)

# Δx, Δy, Δz where Adj non-zero:
ΔX_nz = similar(nonzeros(matrix)); ΔY_nz = similar(nonzeros(matrix)); ΔZ_nz = similar(nonzeros(matrix))


Tmp_D3_nz = similar(nonzeros(matrix))

# 1/(4π * ϵ0 * ϵ_r) * Δx/distance^3, same for Δy and Δz:
S_X = similar(matrix); S_Y = similar(matrix); S_Z = similar(matrix)

# E-field components:
Field_X = similar(X); Field_Y = similar(Y); Field_Z = similar(Z)

## to be put in _drift_charge! after istep for loop
_set_to_zero_vector!(step_vectors)
_add_fieldvector_drift!(step_vectors, current_pos, done, electric_field, det, S)
_add_fieldvector_selfrepulsion!(step_vectors, XI, YI, ZI, XJ, YJ, ZJ, view_XI, view_YI, view_ZI, view_XJ, view_YJ, view_ZJ, ΔX_nz, ΔY_nz, ΔZ_nz, charges, S_X, S_Y, S_Z, Field_X, Field_Y, Field_Z, ϵ_r)