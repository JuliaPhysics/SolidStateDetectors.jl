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

_common_length(dp::EHDriftPath{T} where {T <: SSDFloat})::Int = 
    max(length(dp.timestamps_e), length(dp.timestamps_h))
_common_length(dps::Vector{EHDriftPath{T}} where {T <: SSDFloat})::Int = 
    maximum(_common_length.(dps))

function _common_time(dp::EHDriftPath{T, TT})::TT where {T <: SSDFloat, TT<:RealQuantity}  
    max(last(dp.timestamps_e), last(dp.timestamps_h))
end
_common_time(dps::Vector{<:EHDriftPath}) = 
maximum(_common_time.(dps))

function _common_timestamps(dp::Union{<:EHDriftPath{T}, Vector{<:EHDriftPath{T}}}, Δt) where {T} 
    range(zero(Δt), step = Δt, stop = typeof(Δt)(_common_time(dp)) + Δt)
end

function get_velocity_vector(interpolation_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, point::CartesianPoint{T})::CartesianVector{T} where {T <: SSDFloat}  
    return CartesianVector{T}(interpolation_field(point.x, point.y, point.z))
end

@inline function get_velocity_vector(interpolated_vectorfield, point::CylindricalPoint{T}) where {T <: SSDFloat}
    return CartesianVector{T}(interpolated_vectorfield(point.r, point.φ, point.z))
end


function _drift_charges(detector::SolidStateDetector{T}, grid::Grid{T, 3}, point_types::PointTypes{T, 3},
                        starting_points::Vector{CartesianPoint{T}}, 
                        velocity_field_e::Interpolations.Extrapolation{<:SVector{3}, 3},
                        velocity_field_h::Interpolations.Extrapolation{<:SVector{3}, 3},
                        Δt::RQ; max_nsteps::Int = 2000, verbose::Bool = true)::Vector{EHDriftPath{T}} where {T <: SSDFloat, RQ <: RealQuantity}

    drift_paths::Vector{EHDriftPath{T}} = Vector{EHDriftPath{T}}(undef, length(starting_points)) 

    dt::T = T(to_internal_units(internal_time_unit, Δt))

    for i in eachindex(starting_points)
        drift_path_e::Vector{CartesianPoint{T}} = zeros(CartesianPoint{T}, max_nsteps )#Vector{CartesianPoint{T}}(undef, max_nsteps)
        drift_path_h::Vector{CartesianPoint{T}} = zeros(CartesianPoint{T}, max_nsteps )#Vector{CartesianPoint{T}}(undef, max_nsteps)
        timestamps_e::Vector{T} = Vector{T}(undef, max_nsteps)
        timestamps_h::Vector{T} = Vector{T}(undef, max_nsteps)
        n_e::Int = _drift_charge!(drift_path_e, timestamps_e, detector, point_types, grid, starting_points[i], dt, velocity_field_e, verbose = verbose)
        n_h::Int = _drift_charge!(drift_path_h, timestamps_h, detector, point_types, grid, starting_points[i], dt, velocity_field_h, verbose = verbose)
        drift_paths[i] = EHDriftPath{T, T}( drift_path_e[1:n_e], drift_path_h[1:n_h], timestamps_e[1:n_e], timestamps_h[1:n_h] )
    end

    return drift_paths 
end

function _drift_charge( detector::SolidStateDetector{T}, grid::Grid{T, 3}, point_types::PointTypes{T, 3},
                       starting_point::CartesianPoint{<:SSDFloat}, 
                       velocity_field_e::Interpolations.Extrapolation{SVector{3, T}, 3},
                       velocity_field_h::Interpolations.Extrapolation{SVector{3, T}, 3},
                       Δt::RealQuantity, max_nsteps::Int = 2000, verbose::Bool = true)::Vector{EHDriftPath{T}} where {T <: SSDFloat}
    return _drift_charges(detector, grid, CartesianPoint{T}.(point_types), [starting_point], velocity_field_e, velocity_field_h, T(Δt.val) * unit(Δt), max_nsteps = max_nsteps, verbose = verbose)
end

@inline _convert_vector(pt::CartesianPoint, ::Val{:cylindrical}) = CylindricalPoint(pt)
@inline _convert_vector(pt::CartesianPoint, ::Val{:cartesian}) = pt

function modulate_surface_drift(p::CartesianVector{T})::CartesianVector{T} where {T <: SSDFloat}
    return p
end

function modulate_driftvector(sv::CartesianVector{T}, cp::CartesianPoint{T}, vdv::Vector{AbstractVirtualVolume{T}})::CartesianVector{T} where {T <: SSDFloat}
    for i in eachindex(vdv)
        if in(cp, vdv[i])
            return modulate_driftvector(sv, cp, vdv[i]) 
        end
    end
    return sv
end
"""
    _drift_charge(...)

Before calling this function one should check that `startpos` is inside `det`: `in(startpos, det`
"""
function _drift_charge!(
                            drift_path::Vector{CartesianPoint{T}},
                            timestamps::Vector{T},
                            det::SolidStateDetector{T, S},
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
            last_real_step_index += 1
            current_pos::CartesianPoint{T} = drift_path[istep - 1]
            stepvector::CartesianVector{T} = get_velocity_vector(velocity_field, _convert_vector(current_pos, Val(S))) * Δt
            stepvector = modulate_driftvector(stepvector, current_pos, det.virtual_drift_volumes)
            if geom_round.(stepvector) == null_step 
                done = true 
            end 
            next_pos::CartesianPoint{T} = current_pos + stepvector
            next_pos_cyl::CylindricalPoint{T} = _convert_vector(next_pos, Val(S))
            if next_pos_cyl in point_types || next_pos_cyl in det
                drift_path[istep] = next_pos
                drifttime += Δt
                timestamps[istep] = drifttime
            else
                crossing_pos::CartesianPoint{T}, cd_point_type::UInt8, boundary_index::Int, surface_normal::CartesianVector{T} = get_crossing_pos(det, grid, current_pos, next_pos)
                if cd_point_type == CD_ELECTRODE
                    drift_path[istep] = crossing_pos
                    drifttime += Δt
                    timestamps[istep] = drifttime
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
                    if i == 1000 && verbose @warn("Handling of charge at floating boundary did not work as intended. Start Position (Cart): $startpos") end
                    drift_path[istep] = next_pos
                    drifttime += Δt * (1 - i * T(0.001))
                    timestamps[istep] = drifttime
                    if geom_round.(next_pos - current_pos) == null_step
                        done = true
                    end
                elseif cd_point_type == CD_BULK
                    if verbose @warn ("Internal error for charge starting at $startpos") end
                    drift_path[istep] = current_pos
                    drifttime += Δt
                    timestamps[istep] = drifttime
                    done = true
                else # elseif cd_point_type == CD_OUTSIDE      
                    if verbose @warn ("Internal error for charge starting at $startpos") end
                    drift_path[istep] = current_pos
                    drifttime += Δt
                    timestamps[istep] = drifttime
                    done = true
                end
            end
        end
    end
    return last_real_step_index
end


function is_real(pt::AbstractCoordinatePoint{T})::Bool where {T <: SSDFloat}
    any(v -> isnan(v) || isinf(v), pt[:] )
end
function is_real(pt::AbstractCoordinateVector{T})::Bool where {T <: SSDFloat}
    any(v -> isnan(v) || isinf(v), pt[:] )
end

# Point types for charge drift: Defined in DetectorGeometries/DetectorGeometries.jl
# const CD_ELECTRODE = 0x00
# const CD_OUTSIDE = 0x01
# const CD_BULK = 0x02
# const CD_FLOATING_BOUNDARY = 0x04 # not 0x03, so that one could use bit operations here...

function get_crossing_pos(  detector::SolidStateDetector{T, S}, grid::Grid{T, 3}, point_in::CartesianPoint{T}, point_out::CartesianPoint{T}; 
                            max_n_iter::Int = 500)::Tuple{CartesianPoint{T}, UInt8, Int, CartesianVector{T}} where {T <: SSDFloat, S}
    point_mid::CartesianPoint{T} = T(0.5) * (point_in + point_out)
    cd_point_type::UInt8, contact_idx::Int, surface_normal::CartesianVector{T} = point_type(detector, grid, _convert_vector(point_mid, Val(S)))
    for i in 1:max_n_iter
        if cd_point_type == CD_BULK
            point_in = point_mid
        elseif cd_point_type == CD_OUTSIDE
            point_out = point_mid
        elseif cd_point_type == CD_ELECTRODE
            break
        else #elseif cd_point_type == CD_FLOATING_BOUNDARY
            break
        end
        point_mid = T(0.5) * (point_in + point_out)
        cd_point_type, contact_idx, surface_normal = point_type(detector, grid, _convert_vector(point_mid, Val(S)))
    end
    return point_mid, cd_point_type, contact_idx, surface_normal
end


include("plot_recipes.jl")