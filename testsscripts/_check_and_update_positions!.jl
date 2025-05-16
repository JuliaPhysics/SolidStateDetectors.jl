#using Plots
using SolidStateDetectors: _check_and_update_position!, get_crossing_pos, project_to_plane, modulate_surface_drift, ConstructiveSolidGeometry
using Unitful
using Interpolations
import LegendHDF5IO
import SolidStateDetectors: ConstructiveSolidGeometry

const CD_ELECTRODE = 0x00
const CD_OUTSIDE = 0x01
const CD_BULK = 0x02
const CD_FLOATING_BOUNDARY = 0x04

include("until_istep.jl")

(; step_vectors, current_pos, done, normal, electric_field, S, charges, ϵ_r, CC, drift_path, timestamps, istep, det, grid, point_types, startpos, Δ_t, verbose) = g_state
#function: _check_and_update_position!(step_vectors, current_pos, done, normal, drift_path, timestamps, istep, det, grid, point_types, startpos, Δt, verbose)
istep = 2 # as dummy
g = Grid;

####exlude this loop
for n in eachindex(normal)
    # when the current_pos doesnt change, then save true in done at the index
    done[n] = current_pos[n] == current_pos[n] + step_vectors[n]
    #checks if next point is in the detector
    normal[n] = done[n] || _is_next_point_in_det(current_pos[n] + step_vectors[n], det, point_types)
end

#exclude first if 
if all(normal)
    #all charges are either finished or still inside the detector => drift normally
    current_pos .+= step_vectors
    drift_path[:, istep] .= current_pos
else
    first(findall(.!normal))
    #all charges that would not be inside after the drift step
    for n in findall(.!normal) # every index where false was saved
        crossing_pos::CartesianPoint{T}, cd_point_type::UInt8, surface_normal::CartesianVector{T} =
            get_crossing_pos(det, point_types, copy(current_pos[n]), current_pos[n] + step_vectors[n])
        if cd_point_type == CD_ELECTRODE # if Electron/ Hole crossed Electrode set done true and save the crossing_position  in the drift_path
            done[n] = true
            drift_path[n, istep] = crossing_pos
            current_pos[n] = crossing_pos
        elseif cd_point_type == CD_FLOATING_BOUNDARY
            projected_vector::CartesianVector{T} = CartesianVector{T}(project_to_plane(step_vectors[n], surface_normal))
            projected_vector = modulate_surface_drift(projected_vector)
            next_pos::CartesianPoint{T} = current_pos[n] + projected_vector
            small_projected_vector = projected_vector * T(0.001)
            i::Int = 0
            while i < 1000 && !(next_pos in det.semiconductor)
                next_pos -= small_projected_vector
                i += 1
            end
            if i == 1000
                if verbose
                    @warn("Handling of charge at floating boundary did not work as intended. Start Position (Cart): $(startpos[n])")
                end
                done[n] = true
                continue
            end
            drift_path[n, istep] = next_pos
            step_vectors *= (1 - i * T(0.001))  # scale down the step_vectors for all other charge clouds
            Δt *= (1 - i * T(0.001))            # scale down Δt for all charge clouds
            done[n] = next_pos == current_pos[n]
            current_pos[n] = next_pos
        else # if cd_point_type == CD_BULK or CD_OUTSIDE
            if verbose
                @warn ("Internal error for charge starting at $(startpos[n])")
            end
            done[n] = true
            drift_path[n, istep] = current_pos[n]
        end
    end
    #drift all other charge clouds normally according to the new Δt_min
    for n in findall(normal)
        current_pos[n] += step_vectors[n]
        drift_path[n, istep] = current_pos[n]
    end
end
timestamps[istep] = timestamps[istep-1] + Δ_t
nothing
# end function _check_and_update_position!


function get_crossing_pos(det::SolidStateDetector{T}, point_types::PointTypes{T,3,S}, pts_in::Vector{CartesianPoint{T}}, pts_out::Vector{CartesianPoint{T}};
    max_n_iter::Int=500)::Tuple{CartesianPoint{T},UInt8,CartesianVector{T}} where {T<:SSDFloat,S}

    # check if the points are already in contacts                    
    if pt_in in det.contacts
        return (pt_in, CD_ELECTRODE, CartesianVector{T}(0, 0, 0))
    end

    direction::CartesianVector{T} = normalize(pt_out - pt_in)
    crossing_pos::Tuple{CartesianPoint{T},UInt8,CartesianVector{T}} = (pt_out, CD_OUTSIDE, CartesianVector{T}(0, 0, 0)) # need undef version for this

    # define a Line between pt_in and pt_out
    line::ConstructiveSolidGeometry.Line{T} = ConstructiveSolidGeometry.Line{T}(pt_in, direction)

    # check if the Line intersects with a surface of a Contact
    for contact in det.contacts
        for surf in ConstructiveSolidGeometry.surfaces(contact.geometry)
            # broadcasting this
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

    # if there is no intersection, check if the next point is in (within the tolerance)
    if (crossing_pos[2] & CD_OUTSIDE > 0) && pt_out in det.contacts
        return (pt_out, CD_ELECTRODE, CartesianVector{T}(0, 0, 0))
    end

    crossing_pos
end
