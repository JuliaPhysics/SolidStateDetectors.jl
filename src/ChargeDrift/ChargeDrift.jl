struct DriftPath{T <: SSDFloat}
    path::Vector{<:AbstractCoordinatePoint{RealQuantity}}
    timestamps::Vector{RealQuantity}
end

struct EHDriftPath{T <: SSDFloat, TT <: RealQuantity}
    e_path::Vector{<:AbstractCoordinatePoint{T}}
    h_path::Vector{<:AbstractCoordinatePoint{T}}
    timestamps_e::Vector{TT}
    timestamps_h::Vector{TT}
end

function _common_time(dp::EHDriftPath{T, TT})::TT where {T <: SSDFloat, TT<:RealQuantity}
    max(last(dp.timestamps_e), last(dp.timestamps_h))
end
_common_time(dps::Vector{<:EHDriftPath}) =
maximum(_common_time.(dps))

function _common_timestamps(dp::Union{<:EHDriftPath{T}, Vector{<:EHDriftPath{T}}}, Δt) where {T}
    range(zero(Δt), step = Δt, stop = typeof(Δt)(_common_time(dp)) + Δt)
end

function get_velocity_vector(interpolation_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, pt::CartesianPoint{T})::CartesianVector{T} where {T <: SSDFloat}
    return CartesianVector{T}(interpolation_field(pt.x, pt.y, pt.z))
end

@inline function get_velocity_vector(interpolated_vectorfield, pt::CylindricalPoint{T}) where {T <: SSDFloat}
    return CartesianVector{T}(interpolated_vectorfield(pt.r, pt.φ, pt.z))
end


function _drift_charges(det::SolidStateDetector{T}, grid::Grid{T, 3}, point_types::PointTypes{T, 3},
                        starting_points::Vector{CartesianPoint{T}},
                        velocity_field_e::Interpolations.Extrapolation{<:SVector{3}, 3},
                        velocity_field_h::Interpolations.Extrapolation{<:SVector{3}, 3},
                        Δt::RQ; max_nsteps::Int = 2000, verbose::Bool = true)::Vector{EHDriftPath{T}} where {T <: SSDFloat, RQ <: RealQuantity}

    drift_paths::Vector{EHDriftPath{T}} = Vector{EHDriftPath{T}}(undef, length(starting_points))

    dt::T = T(to_internal_units(Δt))

    for i in eachindex(starting_points)
        drift_path_e::Vector{CartesianPoint{T}} = zeros(CartesianPoint{T}, max_nsteps )#Vector{CartesianPoint{T}}(undef, max_nsteps)
        drift_path_h::Vector{CartesianPoint{T}} = zeros(CartesianPoint{T}, max_nsteps )#Vector{CartesianPoint{T}}(undef, max_nsteps)
        timestamps_e::Vector{T} = Vector{T}(undef, max_nsteps)
        timestamps_h::Vector{T} = Vector{T}(undef, max_nsteps)
        n_e::Int = _drift_charge!(drift_path_e, timestamps_e, det, point_types, grid, starting_points[i], dt, velocity_field_e, verbose = verbose)
        n_h::Int = _drift_charge!(drift_path_h, timestamps_h, det, point_types, grid, starting_points[i], dt, velocity_field_h, verbose = verbose)
        drift_paths[i] = EHDriftPath{T, T}( drift_path_e[1:n_e], drift_path_h[1:n_h], timestamps_e[1:n_e], timestamps_h[1:n_h] )
    end

    return drift_paths
end

function modulate_surface_drift(p::CartesianVector{T})::CartesianVector{T} where {T <: SSDFloat}
    return p
end

function modulate_driftvector(sv::CartesianVector{T}, pt::CartesianPoint{T}, vdv::Vector{AbstractVirtualVolume{T}})::CartesianVector{T} where {T <: SSDFloat}
    for i in eachindex(vdv)
        if in(pt, vdv[i])
            return modulate_driftvector(sv, pt, vdv[i])
        end
    end
    return sv
end
modulate_driftvector(sv::CartesianVector{T}, pt::CartesianPoint{T}, vdv::Missing) where {T} = sv

@inline function _is_next_point_in_det(pt_car::CartesianPoint{T}, pt_cyl::CylindricalPoint{T}, det::SolidStateDetector{T}, point_types::PointTypes{T, 3, Cylindrical})::Bool where {T <: SSDFloat}
    pt_cyl in point_types || pt_cyl in det
end
@inline function _is_next_point_in_det(pt_car::CartesianPoint{T}, pt_cyl::CylindricalPoint{T}, det::SolidStateDetector{T}, point_types::PointTypes{T, 3, Cartesian})::Bool where {T <: SSDFloat}
    pt_car in point_types || pt_car in det
end

function project_to_plane(v⃗::AbstractArray, n⃗::AbstractArray) #Vector to be projected, #normal vector of plane
    # solve (v⃗+λ*n⃗) ⋅ n⃗ = 0
    # n⃗ = n⃗ ./ norm(n⃗)
    λ = -1 * dot(v⃗, n⃗) / dot(n⃗, n⃗)
    SVector{3,eltype(v⃗)}(v⃗[1] + λ * n⃗[1], v⃗[2] + λ * n⃗[2], v⃗[3] + λ * n⃗[3])
end


# """
#     _drift_charge!(...)
# 
# Before calling this function one should check that `startpos` is inside `det`: `in(startpos, det)`
# """
function _drift_charge!(
                            drift_path::Vector{CartesianPoint{T}},
                            timestamps::Vector{T},
                            det::SolidStateDetector{T},
                            point_types::PointTypes{T, 3, S},
                            grid::Grid{T, 3, S},
                            startpos::CartesianPoint{T},
                            Δt::T,
                            velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3};
                            verbose::Bool = true
                        )::Int where {T <: SSDFloat, S}
    drifttime::T = zero(T)
    done::Bool = false
    drift_path[1] = startpos
    timestamps[1] = zero(T)
    null_step::CartesianVector{T} = CartesianVector{T}(0, 0, 0)
    last_real_step_index::Int = 1
    @inbounds for istep in eachindex(drift_path)[2:end]
        if done == false
            current_pos::CartesianPoint{T} = drift_path[istep - 1]
            stepvector::CartesianVector{T} = get_velocity_vector(velocity_field, _convert_point(current_pos, S)) * Δt
            stepvector = modulate_driftvector(stepvector, current_pos, det.virtual_drift_volumes)
            # if geom_round.(stepvector) == null_step
            if stepvector == null_step
                done = true
            end
            next_pos::CartesianPoint{T} = current_pos + stepvector
            next_pos_cyl::CylindricalPoint{T} = CylindricalPoint(next_pos)
            if _is_next_point_in_det(next_pos, next_pos_cyl, det, point_types)
                drift_path[istep] = next_pos
                drifttime += Δt
                timestamps[istep] = drifttime
                last_real_step_index += 1
            else
                crossing_pos::CartesianPoint{T}, cd_point_type::UInt8, boundary_index::Int, surface_normal::CartesianVector{T} = get_crossing_pos(det, grid, current_pos, next_pos)
                if cd_point_type == CD_ELECTRODE
                    drift_path[istep] = crossing_pos
                    drifttime += Δt
                    timestamps[istep] = drifttime
                    last_real_step_index += 1
                    done = true
                elseif cd_point_type == CD_FLOATING_BOUNDARY
                    projected_vector::CartesianVector{T} = CartesianVector{T}(project_to_plane(stepvector, surface_normal))
                    projected_vector = modulate_surface_drift(projected_vector)
                    next_pos = current_pos + projected_vector
                    # ToDo: We actually need a time array as well to do this properly...
                    small_projected_vector = projected_vector * T(0.001)
                    i::Int = 0
                    while i < 1000 && !(next_pos in det)
                        next_pos -= small_projected_vector
                        i += 1
                    end
                    if i == 1000
                        if verbose @warn("Handling of charge at floating boundary did not work as intended. Start Position (Cart): $startpos") end
                        done = true
                        continue
                    end
                    drift_path[istep] = next_pos
                    drifttime += Δt * (1 - i * T(0.001))
                    timestamps[istep] = drifttime
                    last_real_step_index += 1
                    # if geom_round.(next_pos - current_pos) == null_step
                    if next_pos - current_pos == null_step
                        done = true
                    end
                else # elseif cd_point_type == CD_BULK  -- or -- cd_point_type == CD_OUTSIDE
                    if verbose @warn ("Internal error for charge starting at $startpos") end
                    drift_path[istep] = current_pos
                    drifttime += Δt
                    timestamps[istep] = drifttime
                    last_real_step_index += 1
                    done = true
                end
            end
        end
    end
    return last_real_step_index
end

# Point types for charge drift: Defined in DetectorGeometries/DetectorGeometries.jl
# const CD_ELECTRODE = 0x00
# const CD_OUTSIDE = 0x01
# const CD_BULK = 0x02
# const CD_FLOATING_BOUNDARY = 0x04 # not 0x03, so that one could use bit operations here...

function get_crossing_pos(  det::SolidStateDetector{T}, grid::Grid{T, 3, S}, pt_in::CartesianPoint{T}, pt_out::CartesianPoint{T};
                            max_n_iter::Int = 500)::Tuple{CartesianPoint{T}, UInt8, Int, CartesianVector{T}} where {T <: SSDFloat, S}
    pt_mid::CartesianPoint{T} = T(0.5) * (pt_in + pt_out)
    cd_point_type::UInt8, contact_idx::Int, surface_normal::CartesianVector{T} = point_type(det, grid, _convert_point(pt_mid, S))
    for i in 1:max_n_iter
        if cd_point_type == CD_BULK
            pt_in = pt_mid
        elseif cd_point_type == CD_OUTSIDE
            pt_out = pt_mid
        elseif cd_point_type == CD_ELECTRODE
            break
        else #elseif cd_point_type == CD_FLOATING_BOUNDARY
            break
        end
        pt_mid = T(0.5) * (pt_in + pt_out)
        cd_point_type, contact_idx, surface_normal = point_type(det, grid, _convert_point(pt_mid, S))
    end
    return pt_mid, cd_point_type, contact_idx, surface_normal
end
