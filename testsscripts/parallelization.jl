
using Plots
using SolidStateDetectors: elementary_charge, ϵ0, distance_squared, normalize
using Unitful
using Interpolations
import LegendHDF5IO

const CD_ELECTRODE = 0x00
const CD_OUTSIDE = 0x01
const CD_BULK = 0x02
const CD_FLOATING_BOUNDARY = 0x04

include("until_istep.jl")
(; step_vectors, current_pos, done, normal, electric_field, S, charges, ϵ_r, CC, diffusion_length, drift_path, timestamps, istep, det, grid, point_types, startpos, Δ_t, verbose) = g_state

_set_to_zero_vector!(step_vectors)
# _add_fieldvector_drift!(step_vectors, current_pos, done, electric_field, det, S)
# self_repulsion && _add_fieldvector_selfrepulsion!(step_vectors, current_pos, done, charges, ϵ_r)
# _get_driftvectors!(step_vectors, done, Δ_t, det.semiconductor.charge_drift_model, CC)
# diffusion && _add_fieldvector_diffusion!(step_vectors, done, diffusion_length)
# _modulate_driftvectors!(step_vectors, current_pos, det.virtual_drift_volumes)

# function _add_fieldvector_drift!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, done::Vector{Bool}, electric_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, det::SolidStateDetector{T}, ::Type{S})::Nothing where {T, S}
# steps to do before function
# converting current_pos to S
# calculating velocity_vector
    velocity_vector = Vector{CartesianVector{Float32}}(undef, length(step_vectors))# defining before istep
    _set_to_zero_vector!(velocity_vector)
    velocity_vector .= get_velocity_vector.(Ref(electric_field), _convert_point.(current_pos, S))
    step_vectors .= step_vectors .+ velocity_vector .* (.! done)
    nothing
# end

######### geometry important
# function _add_fieldvector_selfrepulsion!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, done::Vector{Bool}, charges::Vector{T}, ϵ_r::T)::Nothing where {T <: SSDFloat}
    #TO DO: ignore charges that are already collected (not trapped though!)

    # possible start:pairs = [(i,j) for i in eachindex(step_vectors), j in eachindex(step_vectors) if j > i]
    for n in eachindex(step_vectors)
        if done[n] continue end
        for m in eachindex(step_vectors)
            if done[m] continue end
            if m > n
                direction::CartesianVector{T} = current_pos[n] .- current_pos[m]
                if iszero(direction) # if the two charges are at the exact same position
                    continue         # don't let them self-repel each other but treat them as same change
                end                  # if diffusion is simulated, they will very likely be separated in the next step
                tmp::T = elementary_charge * inv(4π * ϵ0 * ϵ_r * max(distance_squared(direction), T(1e-10))) # minimum distance is 10µm 
                step_vectors[n] += charges[m] * tmp * normalize(direction)
                step_vectors[m] -= charges[n] * tmp * normalize(direction)
            end
        end
    end
    nothing
# end

# function _get_driftvectors!(step_vectors::Vector{CartesianVector{T}}, done::Vector{Bool}, Δt::T, cdm::AbstractChargeDriftModel{T}, ::Type{Electron})::Nothing where {T <: SSDFloat}
cdm = det.semiconductor.charge_drift_model
    for n in eachindex(step_vectors)
        if !done[n]
            step_vectors[n] = getVe(SVector{3,T}(step_vectors[n]), cdm) * Δt
        end
    end
    nothing
# end

# function _add_fieldvector_diffusion!(step_vectors::Vector{CartesianVector{T}}, done::Vector{Bool}, length::T = T(0.5e3))::Nothing where {T <: SSDFloat}
# steps to do 
# converting step_vectors to a matrix
# step_vectors_matrix = reduce(vcat, [vec' for vec in step_vectors])
sinθ::T, cosθ::T = sincos(acos(T(2*rand() - 1)))
sinφ::T, cosφ::T = sincos(T(rand()*2π))
#calculating the direction of the diffusion
dif_vector = CartesianVector{T}( diffusion_length * cosφ * sinθ, diffusion_length * sinφ * sinθ, diffusion_length * cosθ )
# adding to step_vectors_matrix if done is false
step_vectors .= step_vectors .+ dif_vector' .* (.! done)
nothing 
# end