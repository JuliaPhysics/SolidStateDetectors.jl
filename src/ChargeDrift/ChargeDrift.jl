struct DriftPath{T <: SSDFloat}
    path::Vector{CartesianPoint{T}}
    timestamps::Vector{T}
end

struct EHDriftPath{T <: SSDFloat}
    e_path::Vector{CartesianPoint{T}}
    h_path::Vector{CartesianPoint{T}}
    timestamps_e::Vector{T}
    timestamps_h::Vector{T}
end

function _common_time(dp::EHDriftPath{T})::T where {T <: SSDFloat}
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
                        starting_points::VectorOfArrays{CartesianPoint{T}}, energies::VectorOfArrays{T},
                        electric_field::Interpolations.Extrapolation{<:SVector{3}, 3},
                        Δt::RQ; max_nsteps::Int = 2000, diffusion::Bool = false, self_repulsion::Bool = false, 
                        end_drift_when_no_field::Bool = true, geometry_check::Bool = false, verbose::Bool = true)::Vector{EHDriftPath{T}} where {T <: SSDFloat, RQ <: RealQuantity}

    drift_paths::Vector{EHDriftPath{T}} = Vector{EHDriftPath{T}}(undef, length(flatview(starting_points)))
    dt::T = T(to_internal_units(Δt))
    
    drift_path_counter::Int = 0
    
    for (i, start_points) in enumerate(starting_points)
        
        n_hits::Int = length(start_points)
        charges::Vector{T} = energies[i] ./ to_internal_units(det.semiconductor.material.E_ionisation)
        
        drift_path_e::Array{CartesianPoint{T}, 2} = Array{CartesianPoint{T}, 2}(undef, n_hits, max_nsteps)
        drift_path_h::Array{CartesianPoint{T}, 2} = Array{CartesianPoint{T}, 2}(undef, n_hits, max_nsteps)
        timestamps_e::Vector{T} = Vector{T}(undef, max_nsteps)
        timestamps_h::Vector{T} = Vector{T}(undef, max_nsteps)
        
        n_e::Int = _drift_charge!( drift_path_e, timestamps_e, det, point_types, grid, start_points, -charges, dt, electric_field, Electron; diffusion, self_repulsion, end_drift_when_no_field, geometry_check, verbose )
        n_h::Int = _drift_charge!( drift_path_h, timestamps_h, det, point_types, grid, start_points,  charges, dt, electric_field, Hole; diffusion, self_repulsion, end_drift_when_no_field, geometry_check, verbose )
        
        for i in eachindex(start_points)
            drift_paths[drift_path_counter + i] = EHDriftPath{T}( drift_path_e[i,1:n_e], drift_path_h[i,1:n_h], timestamps_e[1:n_e], timestamps_h[1:n_h] )
        end
        
        drift_path_counter += n_hits
    end
    
    return drift_paths
end

function modulate_surface_drift(p::CartesianVector{T})::CartesianVector{T} where {T <: SSDFloat}
    return p
end

function modulate_driftvector(sv::CartesianVector{T}, pt::CartesianPoint{T}, vdv::Vector{<:AbstractVirtualVolume{T}})::CartesianVector{T} where {T <: SSDFloat}
    for i in eachindex(vdv)
        if in(pt, vdv[i])
            return modulate_driftvector(sv, pt, vdv[i])
        end
    end
    return sv
end
modulate_driftvector(sv::CartesianVector{T}, pt::CartesianPoint{T}, vdv::Missing) where {T} = sv

@inline function _is_next_point_in_det(pt::AbstractCoordinatePoint{T}, det::SolidStateDetector{T}, point_types::PointTypes{T, 3, S})::Bool where {T <: SSDFloat, S}
    _convert_point(pt, S) in point_types || (pt in det.semiconductor && !(pt in det.contacts))
end

function project_to_plane(v⃗::AbstractArray, n⃗::AbstractArray) #Vector to be projected, #normal vector of plane
    # solve (v⃗+λ*n⃗) ⋅ n⃗ = 0
    # n⃗ = n⃗ ./ norm(n⃗)
    λ = -1 * dot(v⃗, n⃗) / dot(n⃗, n⃗)
    SVector{3,eltype(v⃗)}(v⃗[1] + λ * n⃗[1], v⃗[2] + λ * n⃗[2], v⃗[3] + λ * n⃗[3])
end

function _set_to_zero_vector!(v::Vector{CartesianVector{T}})::Nothing where {T <: SSDFloat}
    for n in eachindex(v)
        v[n] = zero(CartesianVector{T})
    end
    nothing
end

function _add_fieldvector_drift!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, done::Vector{Bool}, 
    electric_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, det::SolidStateDetector{T}, ::Type{S}, end_drift_when_no_field::Bool = true)::Nothing where {T, S}
    for n in eachindex(step_vectors)
       if !done[n]
           step_vectors[n] += get_velocity_vector(electric_field, _convert_point(current_pos[n], S))
           done[n] = end_drift_when_no_field && (step_vectors[n] == CartesianVector{T}(0,0,0))
       end
    end
    nothing
end

function _add_fieldvector_diffusion!(step_vectors::Vector{CartesianVector{T}}, done::Vector{Bool}, length::T = T(0.5e3))::Nothing where {T <: SSDFloat}
    for n in eachindex(step_vectors)
        if done[n] continue end
        sinθ::T, cosθ::T = sincos(acos(T(2*rand() - 1)))
        sinφ::T, cosφ::T = sincos(T(rand()*2π))
        step_vectors[n] += CartesianVector{T}( length * cosφ * sinθ, length * sinφ * sinθ, length * cosθ )
    end
    nothing 
end
function _add_fieldvector_diffusion!(step_vectors::Vector{CartesianVector{T}}, done::Vector{Bool}, Δt::T, cdm::AbstractChargeDriftModel{T}, current_pos::Vector{CartesianPoint{T}}, crystal_temperature::T, ::Type{CC})::Nothing where {T <: SSDFloat, CC <: ChargeCarrier}
    for n in eachindex(step_vectors)
        if done[n] continue end

        mu = calculate_mobility(cdm, current_pos[n], CC) # mu in m^2/V/s
        D = mu*crystal_temperature*kB/elementary_charge # D in m^2/s
        length = sqrt(6*D*Δt)

        sinθ::T, cosθ::T = sincos(acos(T(2*rand() - 1)))
        sinφ::T, cosφ::T = sincos(T(rand()*2π))
        step_vectors[n] += CartesianVector{T}( length * cosφ * sinθ, length * sinφ * sinθ, length * cosθ )
    end
    nothing 
end

function _add_fieldvector_selfrepulsion!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, done::Vector{Bool}, charges::Vector{T}, ϵ_r::T)::Nothing where {T <: SSDFloat}
    #TO DO: ignore charges that are already collected (not trapped though!)
    for n in eachindex(step_vectors)
        if done[n] continue end
        for m in eachindex(step_vectors)
            if done[m] continue end
            if m > n
                direction::CartesianVector{T} = current_pos[n] - current_pos[m]
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
end

function _modulate_driftvectors!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, vdv::Vector{V})::Nothing where {T <: SSDFloat, V <: AbstractVirtualVolume{T}}
    for n in eachindex(step_vectors)
        step_vectors[n] = modulate_driftvector(step_vectors[n], current_pos[n], vdv)
    end
    nothing
end
_modulate_driftvectors!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, ::Missing) where {T <: SSDFloat} = nothing

function _get_driftvectors!(step_vectors::Vector{CartesianVector{T}}, done::Vector{Bool}, Δt::T, cdm::AbstractChargeDriftModel{T}, current_pos::Vector{CartesianPoint{T}}, ::Type{Electron})::Nothing where {T <: SSDFloat}
    for n in eachindex(step_vectors)
        if !done[n]
            step_vectors[n] = getVe(SVector{3,T}(step_vectors[n]), cdm, current_pos[n]) * Δt
        end
    end
    nothing
end

function _get_driftvectors!(step_vectors::Vector{CartesianVector{T}}, done::Vector{Bool}, Δt::T, cdm::AbstractChargeDriftModel{T}, current_pos::Vector{CartesianPoint{T}}, ::Type{Hole})::Nothing where {T <: SSDFloat}
    for n in eachindex(step_vectors)
        if !done[n]
            step_vectors[n] = getVh(SVector{3,T}(step_vectors[n]), cdm, current_pos[n]) * Δt
        end
    end
    nothing
end

function _check_and_update_position!(
            step_vectors::Vector{CartesianVector{T}}, 
            current_pos::Vector{CartesianPoint{T}},
            done::Vector{Bool},
            normal::Vector{Bool},
            drift_path::Array{CartesianPoint{T},2},
            timestamps::Vector{T},
            istep::Int,
            det::SolidStateDetector{T},
            g::Grid{T, 3, S},
            point_types::PointTypes{T, 3, S},
            startpos::AbstractVector{CartesianPoint{T}},
            Δt::T,
            geometry_check::Bool,
            verbose::Bool
        )::Nothing where {T <: SSDFloat, S}
        
    for n in eachindex(normal)
        done[n] = current_pos[n] == current_pos[n] + step_vectors[n]
        normal[n] = done[n] || _is_next_point_in_det(current_pos[n]+step_vectors[n], det, point_types)
    end
    
    if all(normal)
        #all charges are either finished or still inside the detector => drift normally
        current_pos .+= step_vectors
        drift_path[:,istep] .= current_pos
    else
        #all charges that would not be inside after the drift step
        for n in findall(.!normal)
            crossing_pos::CartesianPoint{T}, cd_point_type::UInt8, surface_normal::CartesianVector{T} = 
                get_crossing_pos(det, point_types, copy(current_pos[n]), current_pos[n] + step_vectors[n])
            if cd_point_type == CD_ELECTRODE
                if !geometry_check || crossing_pos in det.contacts
                    done[n] = true
                    drift_path[n,istep] = crossing_pos
                    current_pos[n] = crossing_pos
                else
                    cd_point_type = CD_FLOATING_BOUNDARY
                end
            end
            if cd_point_type == CD_FLOATING_BOUNDARY
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
                    if verbose @warn("Handling of charge at floating boundary did not work as intended. Start Position (Cart): $(startpos[n])") end
                    done[n] = true
                    continue
                end
                drift_path[n,istep] = next_pos
                step_vectors *= (1 - i * T(0.001))  # scale down the step_vectors for all other charge clouds
                Δt *= (1 - i * T(0.001))            # scale down Δt for all charge clouds
                done[n] = next_pos == current_pos[n]
                current_pos[n] = next_pos
            elseif cd_point_type!= CD_ELECTRODE # if cd_point_type == CD_BULK or CD_OUTSIDE
                if verbose @warn ("Internal error for charge starting at $(startpos[n])") end
                done[n] = true
                drift_path[n,istep] = current_pos[n]
            end  
        end    
        #drift all other charge clouds normally according to the new Δt_min
        for n in findall(normal)
            current_pos[n] += step_vectors[n]
            drift_path[n,istep] = current_pos[n]
        end
    end
    timestamps[istep] = timestamps[istep-1] + Δt
    nothing
end


function _drift_charge!(
                            drift_path::Array{CartesianPoint{T},2},
                            timestamps::Vector{T},
                            det::SolidStateDetector{T},
                            point_types::PointTypes{T, 3, S},
                            grid::Grid{T, 3, S},
                            startpos::AbstractVector{CartesianPoint{T}},
                            charges::Vector{T},
                            Δt::T,
                            electric_field::Interpolations.Extrapolation{<:StaticVector{3}, 3},
                            ::Type{CC};
                            diffusion::Bool = false,
                            self_repulsion::Bool = false,
                            end_drift_when_no_field::Bool = true,
                            geometry_check::Bool = false,
                            verbose::Bool = true
                        )::Int where {T <: SSDFloat, S, CC <: ChargeCarrier}
                        
    n_hits::Int, max_nsteps::Int = size(drift_path)
    drift_path[:,1] = startpos
    timestamps[1] = zero(T)
    ϵ_r::T = T(det.semiconductor.material.ϵ_r)

    diffusion_length::T = if diffusion
        if CC == Electron && haskey(det.semiconductor.material, :De)
            sqrt(6*_parse_value(T, det.semiconductor.material.De, u"m^2/s") * Δt)
        elseif CC == Hole && haskey(det.semiconductor.material, :Dh)
            sqrt(6*_parse_value(T, det.semiconductor.material.Dh, u"m^2/s") * Δt)
        else 
            @warn "Since v0.9.0, diffusion is modelled via diffusion coefficients `De` (for electrons) and `Dh` (for holes).\n" *
                  "Please update your material properties and pass the diffusion coefficients as `De` and `Dh`.\n" *
                  "You can update it in src/MaterialProperties/MaterialProperties.jl or by overwriting\n" *
                  "`SolidStateDetectors.material_properties` in your julia session and reloading the simulation, e.g.\n
                   SolidStateDetectors.material_properties[:HPGe] = (
                      E_ionisation = 2.95u\"eV\",
                      f_fano = 0.129,
                      ϵ_r = 16.0,
                      ρ = 5.323u\"g*cm^-3\",
                      name = \"High Purity Germanium\",
                      ml = 1.64,
                      mt = 0.0819,
                      De = 200u\"cm^2/s\", # new value 200cm^2/s 
                      Dh = 200u\"cm^2/s\"  # new value 200cm^2/s
                   )\n\n" *
                  "More information can be found at:\n" *
                  "https://juliaphysics.github.io/SolidStateDetectors.jl/stable/man/charge_drift/#Diffusion \n"
            @info "Ignoring diffusion for now"
            diffusion = false
            zero(T)
        end
    else
        zero(T)
    end
    crystal_temperature::T = det.semiconductor.temperature

    last_real_step_index::Int = 1
    current_pos::Vector{CartesianPoint{T}} = deepcopy(startpos)
    step_vectors::Vector{CartesianVector{T}} = Vector{CartesianVector{T}}(undef, n_hits)
    done::Vector{Bool} = broadcast(pt -> !_is_next_point_in_det(pt, det, point_types), startpos)
    normal::Vector{Bool} = deepcopy(done)

    # if the drifting model has a method for calculate_mobility function (scalar mobility), then we use the mobility-tied diffusion model
    cdm = det.semiconductor.charge_drift_model
    ## TODO: make this a keyword argument, so the user can choose if they would like to use this ?
    use_mobility_tied_diffusion = hasmethod(calculate_mobility, Tuple{typeof(cdm), CartesianPoint{T}, Type{CC}})
    
    @inbounds for istep in 2:max_nsteps
        last_real_step_index += 1
        _set_to_zero_vector!(step_vectors)
        _add_fieldvector_drift!(step_vectors, current_pos, done, electric_field, det, S, end_drift_when_no_field)
        self_repulsion && _add_fieldvector_selfrepulsion!(step_vectors, current_pos, done, charges, ϵ_r)
        _get_driftvectors!(step_vectors, done, Δt, cdm, current_pos, CC)
        diffusion && (use_mobility_tied_diffusion ? _add_fieldvector_diffusion!(step_vectors, done, Δt, cdm, current_pos, crystal_temperature, CC) : _add_fieldvector_diffusion!(step_vectors, done, diffusion_length))
        _modulate_driftvectors!(step_vectors, current_pos, det.virtual_drift_volumes)
        _check_and_update_position!(step_vectors, current_pos, done, normal, drift_path, timestamps, istep, det, grid, point_types, startpos, Δt, geometry_check, verbose)
        if all(done) break end
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
    if pt_in in det.contacts return (pt_in, CD_ELECTRODE, zero(CartesianVector{T})) end  
    
    direction::CartesianVector{T} = normalize(pt_out - pt_in)
    crossing_pos::Tuple{CartesianPoint{T}, UInt8, CartesianVector{T}} = (pt_out, CD_OUTSIDE, zero(CartesianVector{T})) # need undef version for this
    
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
    
    # if there is no intersection, check if the next point is in (within the tolerance)
    if (crossing_pos[2] & CD_OUTSIDE > 0) && pt_out in det.contacts 
        return (pt_out, CD_ELECTRODE, CartesianVector{T}(0,0,0)) 
    end

    crossing_pos
end
