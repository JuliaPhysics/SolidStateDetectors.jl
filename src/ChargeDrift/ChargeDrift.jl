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

@inline function get_velocity_vector(velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, pt::CartesianPoint{T})::CartesianVector{T} where {T <: SSDFloat}
    return CartesianVector{T}(velocity_field(pt.x, pt.y, pt.z))
end

@inline function get_velocity_vector(velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, pt::CylindricalPoint{T}) where {T <: SSDFloat}
    return CartesianVector{T}(velocity_field(pt.r, pt.φ, pt.z))
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

@inline function _is_next_point_in_det(pt::CartesianPoint{T}, det::SolidStateDetector{T}, point_types::PointTypes{T, 3, Cylindrical})::Bool where {T <: SSDFloat}
    pt_cyl::CylindricalPoint{T} = CylindricalPoint(pt)
    pt_cyl in point_types || pt_cyl in det.semiconductor
end
@inline function _is_next_point_in_det(pt::CartesianPoint{T}, det::SolidStateDetector{T}, point_types::PointTypes{T, 3, Cartesian})::Bool where {T <: SSDFloat}
    pt in point_types || pt in det.semiconductor
end

function project_to_plane(v⃗::AbstractArray, n⃗::AbstractArray) #Vector to be projected, #normal vector of plane
    # solve (v⃗+λ*n⃗) ⋅ n⃗ = 0
    # n⃗ = n⃗ ./ norm(n⃗)
    λ = -1 * dot(v⃗, n⃗) / dot(n⃗, n⃗)
    SVector{3,eltype(v⃗)}(v⃗[1] + λ * n⃗[1], v⃗[2] + λ * n⃗[2], v⃗[3] + λ * n⃗[3])
end


function _get_stepvector_drift(current_pos::CartesianPoint{T}, ::Type{S}, det::SolidStateDetector{T}, 
                               velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, Δt::T) where {T, S}
    stepvector::CartesianVector{T} = get_velocity_vector(velocity_field, _convert_point(current_pos, S)) * Δt
    stepvector = modulate_driftvector(stepvector, current_pos, det.virtual_drift_volumes)
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
    
    current_pos::CartesianPoint{T} = CartesianPoint{T}(0, 0, 0)
    step_vector::CartesianVector{T} = CartesianVector{T}(0, 0, 0)
    next_pos::CartesianPoint{T} = CartesianPoint{T}(0, 0, 0)
    
    @inbounds for istep in eachindex(drift_path)[2:end]
        if !done
            current_pos = drift_path[istep - 1]
            stepvector = _get_stepvector_drift(current_pos, S, det, velocity_field, Δt)
            done = (stepvector == null_step)
            next_pos = current_pos + stepvector
            
            if _is_next_point_in_det(next_pos, det, point_types)
                drift_path[istep] = next_pos
                drifttime += Δt
                timestamps[istep] = drifttime
                last_real_step_index += 1
                done |= next_pos in det.contacts # end the drift if step ended in a contact
            else
                crossing_pos::CartesianPoint{T}, cd_point_type::UInt8, surface_normal::CartesianVector{T} = get_crossing_pos(det, point_types, current_pos, next_pos)
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
                    done |= (next_pos - current_pos == null_step)
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

# Point types for charge drift
const CD_ELECTRODE = 0x00
const CD_OUTSIDE = 0x01
const CD_BULK = 0x02
const CD_FLOATING_BOUNDARY = 0x04 # not 0x03, so that one could use bit operations here...


function get_crossing_pos(  det::SolidStateDetector{T}, point_types::PointTypes{T, 3, S}, pt_in::CartesianPoint{T}, pt_out::CartesianPoint{T};
                            max_n_iter::Int = 500)::Tuple{CartesianPoint{T}, UInt8, CartesianVector{T}} where {T <: SSDFloat, S}
    
    # check if the points are already in contacts                    
    if pt_in in det.contacts return (pt_in, CD_ELECTRODE, CartesianVector{T}(0,0,0)) end 
    if pt_out in det.contacts return (pt_out, CD_ELECTRODE, CartesianVector{T}(0,0,0)) end 
    
    
    direction::CartesianVector{T} = normalize(pt_out - pt_in)
    crossing_pos::Tuple{CartesianPoint{T}, UInt8, CartesianVector{T}} = (pt_out, CD_OUTSIDE, CartesianVector{T}(0,0,0)) # need undef version for this
    
    # define a Line between pt_in and pt_out
    line::ConstructiveSolidGeometry.Line{T} = ConstructiveSolidGeometry.Line{T}(pt_in, direction)
    
    # check if the Line intersects with a surface of a Contact
    for contact in det.contacts
        for surf in ConstructiveSolidGeometry.surfaces(contact.geometry)
            for pt in ConstructiveSolidGeometry.intersection(surf, line)
                if pt in contact && 
                    0 ≤ (pt - pt_in) ⋅ direction ≤ 1 &&              # pt within pt_in and pt_out
                    norm(pt - pt_in) < norm(crossing_pos[1] - pt_in) # pt closer to pt_in that previous crossing_pos
                    crossing_pos = (pt, CD_ELECTRODE, normalize(ConstructiveSolidGeometry.normal(surf, pt)))
                end
            end
        end
    end
    
    # if the Line does not intersect with a surface of a Contact, check if it does intersect with the surface of the Semiconductor
    if crossing_pos[2] & CD_OUTSIDE > 0 
        tol::T = 5000 * ConstructiveSolidGeometry.csg_default_tol(T)
        # check if the Line intersects with a surface of the Semiconductor
        for surf in ConstructiveSolidGeometry.surfaces(det.semiconductor.geometry)
            for pt in ConstructiveSolidGeometry.intersection(surf, line)
                normal = normalize(ConstructiveSolidGeometry.normal(surf, pt))
                if pt + tol * normal in det.semiconductor &&         # point "before" crossing_pos should be in
                   !(pt - tol * normal in det.semiconductor) &&      # point "after" crossing_pos should be out
                    0 ≤ (pt - pt_in) ⋅ direction ≤ 1 &&              # pt within pt_in and pt_out
                    norm(pt - pt_in) < norm(crossing_pos[1] - pt_in) # pt closer to pt_in that previous crossing_pos
                    CD_POINTTYPE::UInt8 = point_types[pt] & update_bit == 0 ? CD_ELECTRODE : CD_FLOATING_BOUNDARY 
                    crossing_pos = (pt, CD_POINTTYPE, normal)
                end
            end 
        end
    end

    crossing_pos
end
