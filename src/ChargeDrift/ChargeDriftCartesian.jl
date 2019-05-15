
"""
    Before calling this function one should check that `startpos` is inside `det`: `in(startpos, det`

"""
function drift_charge!(
                            drift_path::Vector{CartesianPoint{T}},
                            det::SolidStateDetector{T, :cartesian},
                            point_types::PointTypes{T, 3, :cartesian},
                            grid::Grid{T, 3, :cartesian},
                            startpos::CartesianPoint{T},
                            Δt::T,
                            velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3};
                            verbose::Bool = true
                        )::Nothing where {T <: SSDFloat}
    done::Bool = false
    drift_path[1] = startpos
    null_step::CartesianVector{T} = CartesianVector{T}(0, 0, 0)
    @inbounds for istep in eachindex(drift_path)[2:end] #end] 
        if done == false
            current_pos::CartesianPoint{T} = drift_path[istep - 1]
            stepvector::CartesianVector{T} = get_velocity_vector(velocity_field, current_pos) * Δt
            if stepvector == null_step done = true end 
            next_pos::CartesianPoint{T} = current_pos + stepvector
            if next_pos in point_types
                drift_path[istep] = next_pos
            elseif next_pos in det
                drift_path[istep] = next_pos
            else
                crossing_pos::CartesianPoint{T}, cd_point_type::UInt8, boundary_index::Int, surface_normal::CartesianVector{T} = get_crossing_pos(det, grid, current_pos, next_pos)
                if cd_point_type == CD_ELECTRODE
                    done = true
                    drift_path[istep] = crossing_pos
                elseif cd_point_type == CD_FLOATING_BOUNDARY
                    projected_vector = project_to_plane(stepvector, surface_normal )
                    next_pos = current_pos + projected_vector
                    # ToDo: We actually need a time array as well to do this properly...
                    small_projected_vector = projected_vector * T(0.001)
                    i::Int = 0
                    while i < 1000 && (next_pos in point_types || !(next_pos in det)) 
                        next_pos -= small_projected_vector
                        i += 1
                    end
                    if i == 1000 && verbose @warn("Handling of charge at floating boundary did not work as intended. Start Position (Cart): $startpos") end
                    drift_path[istep] = next_pos
                elseif cd_point_type == CD_BULK
                    if verbose @warn ("Internal error for charge stating at $startpos") end
                    drift_path[istep] = current_pos
                    done = true
                else # elseif cd_point_type == CD_OUTSIDE      
                    if verbose @warn ("Internal error for charge stating at $startpos") end
                    drift_path[istep] = current_pos
                    done = true
                end
            end
        elseif done == true
            drift_path[istep] = drift_path[istep-1]
        end
    end
    return nothing
end

# Point types for charge drift: Defined in DetectorGeometries/DetectorGeometries.jl
# const CD_ELECTRODE = 0x00
# const CD_OUTSIDE = 0x01
# const CD_BULK = 0x02
# const CD_FLOATING_BOUNDARY = 0x04 # not 0x03, so that one could use bit operations here...

function get_crossing_pos(  detector::SolidStateDetector{T, :cartesian}, grid::Grid{T, 3}, point_in::CartesianPoint{T}, point_out::CartesianPoint{T}; 
                            max_n_iter::Int = 2000)::Tuple{CartesianPoint{T}, UInt8, Int, CartesianVector{T}} where {T <: SSDFloat}
    point_mid::CartesianPoint{T} = T(0.5) * (point_in + point_out)
    cd_point_type::UInt8, contact_idx::Int, surface_normal::CartesianVector{T} = point_type(detector, grid, point_mid)
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
        cd_point_type, contact_idx, surface_normal = point_type(detector, grid, point_mid)
    end
    return point_mid, cd_point_type, contact_idx, surface_normal
end

