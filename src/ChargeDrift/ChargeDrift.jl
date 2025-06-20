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

struct EDepEvent{T<:RealQuantity, U<:RealQuantity}
    Pos::AbstractVector{<:CartesianVector{T}}
    Q::AbstractVector{<:U}
end

abstract type ChargeCarrier end
abstract type Electron <: ChargeCarrier end 
abstract type Hole <: ChargeCarrier end

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


function calculate_drift_fields!(sim::Simulation{T};
    use_nthreads::Int = Base.Threads.nthreads())::Nothing where {T <: SSDFloat}
    @warn "Since v0.7.0, drift fields do not need to be calculated anymore.\n`calculate_drift_fields!(sim)` can be removed."
    nothing
end
@deprecate apply_charge_drift_model!(args...; kwargs...) calculate_drift_fields!(args...; kwargs...)

#originally in SolidStateDetectors.jl
function drift_charges( sim::Simulation{T}, starting_positions::VectorOfArrays{CartesianPoint{T}}, energies::VectorOfArrays{T};
    Δt::RealQuantity = 5u"ns", max_nsteps::Int = 1000, diffusion::Bool = false, self_repulsion::Bool = false, verbose::Bool = true )::Vector{EHDriftPath{T}} where {T <: SSDFloat}
    return _drift_charges(   sim.detector, sim.point_types.grid, sim.point_types, starting_positions, energies, 
         interpolated_vectorfield(sim.electric_field), Δt, 
         max_nsteps = max_nsteps, diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose)
end

# for one event
# function drift_charges( sim::Simulation{T}, starting_positions::AbstractVector{CartesianPoint{T}}, energies::AbstractVector{<:Quantity{T}};
#                         Δt::RealQuantity = 5u"ns", max_nsteps::Int = 1000, diffusion::Bool = false, self_repulsion::Bool = false, verbose::Bool = true ) where {T <: SSDFloat} # ::Vector{EHDriftPath{T}}
#     starting_positions = convert_to_cartesian.(starting_positions)
#     starting_positions = [starting_positions]
#     energies = [energies]
#     edep_events = StructVector{EDepEvent}((starting_positions, energies))
#     return drift_charges(sim, edep_events; 
#         Δt = Δt, max_nsteps = max_nsteps, diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose)
# end

function drift_charges( sim::Simulation{T}, starting_positions::AbstractVector{CartesianPoint{T}}, energies::AbstractVector{<:Quantity{T}};
    Δt::RealQuantity = 5u"ns", max_nsteps::Int = 1000, diffusion::Bool = false, self_repulsion::Bool = false, verbose::Bool = true ) where {T <: SSDFloat} #::Vector{EHDriftPath{T}}
    starting_positions = [starting_positions]
    energies = [energies]
    edep_events = StructVector{EDepEvent}((starting_positions, energies))
    return drift_charges(sim, edep_events; 
    Δt = Δt, max_nsteps = max_nsteps, diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose)
end

# for event table
function drift_charges( sim::Simulation{T}, edep_events_table::Any;
    Δt::RealQuantity = 5u"ns", max_nsteps::Int = 1000, diffusion::Bool = false, self_repulsion::Bool = false, verbose::Bool = true )::Vector{EHDriftPath{T}} where {T <: SSDFloat} 
    edep_events = convert_edepevents(edep_events_table)
    return drift_charges(sim, edep_events; Δt = Δt, max_nsteps = max_nsteps, diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose)
end

# strips Unit asnd calls _drift_charges
function drift_charges(sim::Simulation{T}, edep_events::StructVector{EDepEvent};
    Δt::RealQuantity = 5u"ns", max_nsteps::Int = 1000, diffusion::Bool = false, self_repulsion::Bool = false, verbose::Bool = true )::Vector{EHDriftPath{T}} where {T <: SSDFloat}
    starting_positions = VectorOfVectors(_ustrip_recursive(edep_events.Pos))
    energies = VectorOfVectors(_ustrip_recursive(edep_events.Q))
    dt::T = to_internal_units(Δt)
    return _drift_charges(sim.detector, sim.point_types.grid, sim.point_types, starting_positions, energies, interpolated_vectorfield(sim.electric_field), dt, max_nsteps = max_nsteps, diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose)
end

function get_signal(sim::Simulation{T, CS}, drift_paths::Vector{EHDriftPath{T}}, energy_depositions::Vector{T}, contact_id::Int; Δt::TT = T(5) * u"ns") where {T <: SSDFloat, CS, TT}
    dt::T = to_internal_units(Δt)
    wpot::Interpolations.Extrapolation{T, 3} = interpolated_scalarfield(sim.weighting_potentials[contact_id])
    timestamps = _common_timestamps( drift_paths, dt )
    signal::Vector{T} = zeros(T, length(timestamps))
    add_signal!(signal, timestamps, drift_paths, energy_depositions, wpot, CS, sim.detector.semiconductor.charge_trapping_model)
    unitless_energy_to_charge = _convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)
    return RDWaveform( range(zero(T) * unit(Δt), step = T(ustrip(Δt)) * unit(Δt), length = length(signal)), signal * unitless_energy_to_charge)
end

# change 
function _drift_charges(det::SolidStateDetector{T}, grid::Grid{T, 3}, point_types::PointTypes{T, 3},
                        starting_points::VectorOfVectors{CartesianPoint{T}}, energies::VectorOfVectors{T},
                        electric_field::Interpolations.Extrapolation{<:SVector{3}, 3},
                        Δt::RQ; max_nsteps::Int = 2000, diffusion::Bool = false, self_repulsion::Bool = false, verbose::Bool = true)::Vector{EHDriftPath{T}} where {T <: SSDFloat, RQ <: RealQuantity} 

    starting_points_pointer = starting_points.elem_ptr
    start_points = flatview(starting_points)
    drift_paths::Vector{EHDriftPath{T}} = Vector{EHDriftPath{T}}(undef, length(start_points))
    dt::T = T(to_internal_units(Δt))
    
    drift_path_counter::Int = 0
    
        
        n_hits::Int = length(start_points)
        charges::Vector{T} = flatview(energies) ./ to_internal_units(det.semiconductor.material.E_ionisation)
        
        drift_path_e::Array{CartesianPoint{T}, 2} = Array{CartesianPoint{T}, 2}(undef, n_hits, max_nsteps)
        drift_path_h::Array{CartesianPoint{T}, 2} = Array{CartesianPoint{T}, 2}(undef, n_hits, max_nsteps)
        timestamps_e::Vector{T} = Vector{T}(undef, max_nsteps)
        timestamps_h::Vector{T} = Vector{T}(undef, max_nsteps)

        n_e = _drift_charge!( drift_path_e, timestamps_e, det, point_types, grid, start_points, -charges, dt, electric_field, Electron, diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose ) # ::Int
        n_h = _drift_charge!( drift_path_h, timestamps_h, det, point_types, grid, start_points,  charges, dt, electric_field, Hole, diffusion = diffusion, self_repulsion = self_repulsion, verbose = verbose ) # ::Int
    
        for i in eachindex(start_points)
            drift_paths[drift_path_counter + i] = EHDriftPath{T}( drift_path_e[i,1:n_e], drift_path_h[i,1:n_h], timestamps_e[1:n_e], timestamps_h[1:n_h] )
        end
        
        drift_path_counter += n_hits

    return drift_paths
end
# region other functions 
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
    #new for flatvector
    @inline function _is_next_point_in_det(pt::CartesianVector{T}, det::SolidStateDetector{T}, point_types::PointTypes{T, 3, S})::Bool where {T <: SSDFloat, S}
        _convert_point(pt, S) in point_types || (pt in det.semiconductor && !(pt in det.contacts))
    end

    function project_to_plane(v⃗::AbstractArray, n⃗::AbstractArray) #Vector to be projected, #normal vector of plane
        # solve (v⃗+λ*n⃗) ⋅ n⃗ = 0
        # n⃗ = n⃗ ./ norm(n⃗)
        λ = -1 * dot(v⃗, n⃗) / dot(n⃗, n⃗)
        SVector{3,eltype(v⃗)}(v⃗[1] + λ * n⃗[1], v⃗[2] + λ * n⃗[2], v⃗[3] + λ * n⃗[3])
    end
        
    function get_cutoff(all_positions, el_field_itp::Interpolations.Extrapolation{<:StaticVector{3}, 3})
        T = eltype(eltype(all_positions))
        arr = Vector{typeof(zero(T) * u"m")}()
        N = length(all_positions)
        k_e = 1 / (4 * π * ϵ0*u"F/m")
        for i in 1:N, j in 1:N
            if j != i
                q1 = elementary_charge * u"C"
                q2 = elementary_charge * u"C"
                efield = el_field_itp(all_positions[i].x, all_positions[i].y, all_positions[i].z) * u"V/m"
                F_el = q1 * efield 
                
                x1 = all_positions[i].x * u"mm"
                y1 = all_positions[i].y * u"mm"
                z1 = all_positions[i].z * u"mm"

                x2 = all_positions[j].x * u"mm"
                y2 = all_positions[j].y * u"mm"
                z2 = all_positions[j].z * u"mm"

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
        return ustrip(isempty(arr) ? zero(eltype(arr)) : maximum(arr))
    end

    function build_cutoff_matrix(current_pos, electric_field::Interpolations.Extrapolation{<:StaticVector{3, T}, 3}) where {T <: SSDFloat}
        d_cutoff = get_cutoff(current_pos, electric_field)
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

    function _set_to_zero_vector!(v::Vector{CartesianVector{T}})::Nothing where {T <: SSDFloat}
        v .= (CartesianVector{T}(0,0,0),)
        # for n in eachindex(v)
        #     v[n] = CartesianVector{T}(0,0,0)
        # end
        nothing
    end

    function _add_fieldvector_drift!(step_vectors::Vector{CartesianVector{T}}, velocity_vector::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, done::Vector{Bool}, electric_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, det::SolidStateDetector{T}, ::Type{S})::Nothing where {T, S}
        _set_to_zero_vector!(velocity_vector)
        velocity_vector .= get_velocity_vector.(Ref(electric_field), _convert_point.(current_pos, S))
        step_vectors .= step_vectors .+ velocity_vector .* (.! done)
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
            nothing
        end
    end

    function _add_fieldvector_selfrepulsion!(step_vectors::AbstractVector{CartesianVector{T}}, 
        XI::AbstractVector{T}, YI::AbstractVector{T}, ZI::AbstractVector{T}, 
        XJ::AbstractVector{T}, YJ::AbstractVector{T}, ZJ::AbstractVector{T}, 
        view_XI::AbstractVector{T}, view_YI::AbstractVector{T}, view_ZI::AbstractVector{T}, 
        view_XJ::AbstractVector{T}, view_YJ::AbstractVector{T}, view_ZJ::AbstractVector{T}, 
        ΔX_nz::AbstractVector{T}, ΔY_nz::AbstractVector{T}, ΔZ_nz::AbstractVector{T}, 
        Tmp_D3_nz, 
        charges::AbstractVector{T}, 
        S_X, S_Y, S_Z, 
        Field_X::AbstractVector{T}, Field_Y::AbstractVector{T}, Field_Z::AbstractVector{T}, 
        ϵ_r::T)::Nothing where {T <: SSDFloat}

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


    function _modulate_driftvectors!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, vdv::Vector{V})::Nothing where {T <: SSDFloat, V <: AbstractVirtualVolume{T}}
        for n in eachindex(step_vectors)
            step_vectors[n] = modulate_driftvector(step_vectors[n], current_pos[n], vdv)
        end
        nothing
    end
    _modulate_driftvectors!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, ::Missing) where {T <: SSDFloat} = nothing

    function _get_driftvectors!(step_vectors::Vector{CartesianVector{T}}, done::Vector{Bool}, Δt::T, cdm::AbstractChargeDriftModel{T}, ::Type{Electron})::Nothing where {T <: SSDFloat}
        for n in eachindex(step_vectors)
            if !done[n]
                step_vectors[n] = getVe(SVector{3,T}(step_vectors[n]), cdm) * Δt
            end
        end
        nothing
    end

    function _get_driftvectors!(step_vectors::Vector{CartesianVector{T}}, done::Vector{Bool}, Δt::T, cdm::AbstractChargeDriftModel{T}, ::Type{Hole})::Nothing where {T <: SSDFloat}
        for n in eachindex(step_vectors)
            if !done[n]
                step_vectors[n] = getVh(SVector{3,T}(step_vectors[n]), cdm) * Δt
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
                    done[n] = true
                    drift_path[n,istep] = crossing_pos
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
                        if verbose @warn("Handling of charge at floating boundary did not work as intended. Start Position (Cart): $(startpos[n])") end
                        done[n] = true
                        continue
                    end
                    drift_path[n,istep] = next_pos
                    step_vectors *= (1 - i * T(0.001))  # scale down the step_vectors for all other charge clouds
                    Δt *= (1 - i * T(0.001))            # scale down Δt for all charge clouds
                    done[n] = next_pos == current_pos[n]
                    current_pos[n] = next_pos
                else # if cd_point_type == CD_BULK or CD_OUTSIDE
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
# endregion

#!!!!!!!!!!!! REMOVE THIS
function _drift_charge_old!(
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

    last_real_step_index::Int = 1
    current_pos::Vector{CartesianPoint{T}} = deepcopy(startpos)
    step_vectors::Vector{CartesianVector{T}} = Vector{CartesianVector{T}}(undef, n_hits)
    done::Vector{Bool} = broadcast(pt -> !_is_next_point_in_det(pt, det, point_types), startpos)
    normal::Vector{Bool} = deepcopy(done)
    
    @inbounds for istep in 2:max_nsteps
        last_real_step_index += 1
        _set_to_zero_vector!(step_vectors)
        _add_fieldvector_drift!(step_vectors, current_pos, done, electric_field, det, S)
        self_repulsion && _add_fieldvector_selfrepulsion!(step_vectors, current_pos, done, charges, ϵ_r)
        _get_driftvectors!(step_vectors, done, Δt, det.semiconductor.charge_drift_model, CC)
        diffusion && _add_fieldvector_diffusion!(step_vectors, done, diffusion_length)
        _modulate_driftvectors!(step_vectors, current_pos, det.virtual_drift_volumes)
        _check_and_update_position!(step_vectors, current_pos, done, normal, drift_path, timestamps, istep, det, grid, point_types, startpos, Δt, verbose)
        if all(done) break end
    end

    return last_real_step_index
end

function _add_fieldvector_drift!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, done::Vector{Bool}, 
    electric_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, det::SolidStateDetector{T}, ::Type{S})::Nothing where {T, S}
    for n in eachindex(step_vectors)
       if !done[n]
           step_vectors[n] += get_velocity_vector(electric_field, _convert_point(current_pos[n], S))
           done[n] = (step_vectors[n] == CartesianVector{T}(0,0,0))
       end
    end
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
                            verbose::Bool = true, 
                        )::Int where {T <: SSDFloat, S, CC <: ChargeCarrier} # ::AbstractVector 
                        
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

    last_real_step_index::Int = 1
    current_pos::Vector{CartesianPoint{T}} = deepcopy(startpos)
    step_vectors::Vector{CartesianVector{T}} = Vector{CartesianVector{T}}(undef, n_hits)
    done::Vector{Bool} = broadcast(pt -> !_is_next_point_in_det(pt, det, point_types), startpos)
    normal::Vector{Bool} = deepcopy(done)

    # Memory allocation for parallelizable _add_fieldvector_selfrepulsion!
    velocity_vector = Vector{CartesianVector{T}}(undef, length(step_vectors))

    matrix, X, Y, Z = build_cutoff_matrix(current_pos, electric_field)

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

        
    @inbounds for istep in 2:max_nsteps
        last_real_step_index += 1
        _set_to_zero_vector!(step_vectors)
        _add_fieldvector_drift!(step_vectors, velocity_vector, current_pos, done, electric_field, det, S)
        self_repulsion && _add_fieldvector_selfrepulsion!(step_vectors, XI, YI, ZI, XJ, YJ, ZJ, view_XI, view_YI, view_ZI, view_XJ, view_YJ, view_ZJ, ΔX_nz, ΔY_nz, ΔZ_nz, Tmp_D3_nz, charges, S_X, S_Y, S_Z, Field_X, Field_Y, Field_Z, ϵ_r)
        _get_driftvectors!(step_vectors, done, Δt, det.semiconductor.charge_drift_model, CC)
        diffusion && _add_fieldvector_diffusion!(step_vectors, done, diffusion_length)
        _modulate_driftvectors!(step_vectors, current_pos, det.virtual_drift_volumes)
        _check_and_update_position!(step_vectors, current_pos, done, normal, drift_path, timestamps, istep, det, grid, point_types, startpos, Δt, verbose)
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
    if pt_in in det.contacts return (pt_in, CD_ELECTRODE, CartesianVector{T}(0,0,0)) end  
    
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
    
    # if there is no intersection, check if the next point is in (within the tolerance)
    if (crossing_pos[2] & CD_OUTSIDE > 0) && pt_out in det.contacts 
        return (pt_out, CD_ELECTRODE, CartesianVector{T}(0,0,0)) 
    end

    crossing_pos
end
