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


#!!!!!!!!!!!!!!!!! Remove either get_velocity_vector or get_velocity_vector_withconv
@inline function get_velocity_vector(velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, pt::CartesianPoint{T})::CartesianVector{T} where {T <: SSDFloat}
    return CartesianVector{T}(velocity_field(pt.x, pt.y, pt.z))
end

@inline function get_velocity_vector(velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, pt::CylindricalPoint{T}) where {T <: SSDFloat}
    return CartesianVector{T}(velocity_field(pt.r, pt.φ, pt.z))
end

@inline function get_velocity_vector_withconv(velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, pt::AbstractCoordinatePoint{T}, ::Type{Cartesian}, done::Bool)::CartesianVector{T} where {T <: SSDFloat}
    if done
        return CartesianVector{T}(0, 0, 0)
    else
        conv_pt = _convert_point(pt, Cartesian)
        return CartesianVector{T}(velocity_field(conv_pt.x, conv_pt.y, conv_pt.z))
    end
end

@inline function get_velocity_vector_withconv(velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, pt::AbstractCoordinatePoint{T}, ::Type{Cylindrical}, done::Bool) where {T <: SSDFloat}
    if done
        return CartesianVector{T}(0, 0, 0)
    else
        conv_pt = _convert_point(pt, Cylindrical)
        return CartesianVector{T}(velocity_field(conv_pt.r, conv_pt.φ, conv_pt.z))
    end
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

    function build_cutoff_matrix(pos, electric_field::Interpolations.Extrapolation{<:StaticVector{3, T}, 3}) where {T <: SSDFloat}
        d_cutoff = get_cutoff(pos, electric_field)
        N = length(pos)

        (;x, y, z) = pos

        row_indices = Int[]
        col_indices = Int[]
        values = Bool[]

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
        return sparse(row_indices, col_indices, fill(one(T), length(row_indices)), N, N)
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

    function _add_fieldvector_selfrepulsion!(field_vectors::AbstractVector{CartesianVector{T}}, current_pos::AbstractVector{CartesianPoint{T}}, done::AbstractVector{Bool}, charges::AbstractVector{T}, ϵ_r::T)::Nothing where {T <: SSDFloat}
        #TO DO: ignore charges that are already collected (not trapped though!)
        for n in eachindex(field_vectors)
            if done[n] continue end
            for m in eachindex(field_vectors)
                if done[m] continue end
                if m > n
                    direction::CartesianVector{T} = current_pos[n] - current_pos[m]
                    if iszero(direction) # if the two charges are at the exact same position
                        continue         # don't let them self-repel each other but treat them as same change
                    end                  # if diffusion is simulated, they will very likely be separated in the next step
                    tmp::T = elementary_charge * inv(4π * ϵ0 * ϵ_r * max(distance_squared(direction), T(1e-10))) # minimum distance is 10µm 
                    field_vectors[n] += charges[m] * tmp * normalize(direction)
                    field_vectors[m] -= charges[n] * tmp * normalize(direction)
                end
            end
            nothing
        end
    end


    function _modulate_driftvectors!(step_vectors::AbstractVector{CartesianVector{T}}, current_pos::AbstractVector{CartesianPoint{T}}, vdv::AbstractVector{V})::Nothing where {T <: SSDFloat, V <: AbstractVirtualVolume{T}}
        for n in eachindex(step_vectors)
            step_vectors[n] = modulate_driftvector(step_vectors[n], current_pos[n], vdv)
        end
        nothing
    end
    _modulate_driftvectors!(step_vectors::AbstractVector{CartesianVector{T}}, current_pos::AbstractVector{CartesianPoint{T}}, ::Missing) where {T <: SSDFloat} = nothing


    function _get_drift_step(field_vector::CartesianVector{T}, cdm::AbstractChargeDriftModel{T}, Δt::T, done::Bool) where {T}
        if !done
            return (getVe(SVector{3,T}(field_vector), cdm) * Δt)::SVector{3,T}
        else
            return SVector{3,T}(0, 0, 0)
        end
    end

    function _get_drift_steps!(
        step_vectors::AbstractVector{<:CartesianVector{T}}, field_vectors::AbstractVector{<:CartesianVector{T}},
        done::AbstractVector{Bool}, Δt::T, cdm::AbstractChargeDriftModel{T}, ::Type{CC}
    )::Nothing where {T <: SSDFloat, CC <: ChargeCarrier}
        broadcast!(_get_drift_step, step_vectors, field_vectors, Ref(cdm), Δt, done)
        nothing
    end


    function _check_and_update_position!(
                step_vectors::AbstractVector{CartesianVector{T}},
                current_pos::AbstractVector{CartesianPoint{T}},
                done::AbstractVector{Bool},
                normal::AbstractVector{Bool},
                drift_path::AbstractArray{CartesianPoint{T},2},
                timestamps::AbstractVector{T},
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


function _get_diffusion_length(det::SolidStateDetector{T}, Δt::Real, ::Type{CC})::T where {T,CC}
    if CC == Electron && haskey(det.semiconductor.material, :De)
        return sqrt(6*_parse_value(T, det.semiconductor.material.De, u"m^2/s") * Δt)
    elseif CC == Hole && haskey(det.semiconductor.material, :Dh)
        return sqrt(6*_parse_value(T, det.semiconductor.material.Dh, u"m^2/s") * Δt)
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
        return zero(T)
    end
end


struct ChargeInteractionState{
    T<:Real, TV <: AbstractVector{T}, TVView <: AbstractVector{T}, TM <: AbstractMatrix{T},
    PV <: StructVector{<:CartesianPoint{T}}, PVView <: StructVector{<:CartesianPoint{T}},
    VV <: StructVector{<:CartesianVector{T}},
    TS <: NamedTuple{(:x,:y,:z), <:Tuple{TM,TM,TM}}
}
    M_adj::TM # sparse charge interaction adjacency matrix, elements must be one or zero, diagonal must be zero
    adj_row_sums::TV # count of non-zero elements in each row of M_adj
    pos::PV # positions of charges
    pos_I::PV # pos expanded to interaction-target adj. indices
    pos_J::PV # pos expanded to interaction-source adj. indices
    pos_vI::PVView # pos_J as a view into the original charge positions
    pos_vJ::PVView # pos_I as a view into the original charge positions
    Δpos_IJ::VV # pos_I - pos_J
    charges::TV # charge values
    charges_J::TV # charges expanded to interaction-source adj. indices
    charges_vJ::TVView # charges_J as a view into the original charges
    Tmp_D3_nz::TV
    S::TS # 1/(4π * ϵ0 * ϵ_r) * [Δx,Δy,Δz]/distance^3
    ϵ_r::T # Relative permittivity
    onesT::TV # Vector of ones, same length as charges
    contr_thresh::TV # Interaction threshold by target charge
    contr_thresh_I::TV # Interaction threshold by target charge, expanded to interaction-target adj. indices
    contr_thresh_vI::TVView # Threshold for interaction for each row of M_adj
end

function ChargeInteractionState(
    pos::StructVector{<:CartesianPoint{T}}, charges::AbstractVector{T}, ϵ_r::T,
    M_adj::AbstractMatrix{T}
) where T
    # ToDo:
    # For `idxJ, idxJ, nzA = findnz(A)`, `nzA` is a copy, while `nonzeros(A)` returns
    # the inner non-zero data of A. In our charge-interaction implementation we depend
    # on them being euqal (equivalent in respect to `idxJ, idxJ`). This seems to be the
    # case for `SparseMatrixCSC` and `CuSparseMatrixCSC`, but is it guaranteed in
    # general? If not, we should test this in this constructor, for the given sparse
    # matrix type.

    n_charges = length(pos)

    issparse(M_adj) || throw(ArgumentError("M_adj must be a sparse matrix"))
    iszero(diag(M_adj)) || throw(ArgumentError("Diagonal of M_adj must be zero"))

    adj_I, adj_J = findnz(M_adj)
    adj_nz = nonzeros(M_adj)
    adj_nnz = length(adj_nz)
    @assert length(adj_I) == length(adj_J) == adj_nnz

    onesT = one.(charges)
    adj_row_sums = M_adj * onesT

    pos_I = similar(pos, adj_nnz)
    @inbounds pos_vI = view(pos, adj_I)
    pos_J = similar(pos, adj_nnz)
    @inbounds pos_vJ = view(pos, adj_J)
    Δpos_IJ = StructVector{CartesianVector{T}}(x = similar(pos_I.x), y = similar(pos_I.y), z = similar(pos_I.z))

    charges_J = similar(charges, adj_nnz)
    @inbounds charges_vJ = view(charges, adj_J)

    Tmp_D3_nz = similar(adj_nz)
    S = (x = _linked_sparsematrix(M_adj), y = _linked_sparsematrix(M_adj), z = _linked_sparsematrix(M_adj))

    contr_thresh = similar(charges, n_charges)
    contr_thresh_I = similar(contr_thresh, adj_nnz)
    @inbounds contr_thresh_vI = view(contr_thresh, adj_I)

    return ChargeInteractionState(
        M_adj, adj_row_sums, pos, pos_I, pos_J, pos_vI, pos_vJ, Δpos_IJ,
        charges, charges_J, charges_vJ,
        Tmp_D3_nz, S, ϵ_r, onesT, contr_thresh, contr_thresh_I, contr_thresh_vI
    )
end

function _linked_sparsematrix(M::SparseMatrixCSC)
    m, n = M.m, M.n
    new_nzval = deepcopy(M.nzval)
    return SparseMatrixCSC(m, n, M.colptr, M.rowval, new_nzval)
end

function _linked_resize!!(C::SparseMatrixCSC, M::SparseMatrixCSC)
    m, n = M.m, M.n
    new_nzval = resize!(C.nzval, length(M.nzval))
    return SparseMatrixCSC(m, n, M.colptr, M.rowval, new_nzval)
end


function Adapt.adapt_structure(adapter, ci_state::ChargeInteractionState{T}) where T <: Real
    M_adj = adapt(adapter, ci_state.M_adj)

    adj_I, adj_J = findnz(M_adj)
    adj_nz = nonzeros(M_adj)
    adj_nnz = length(adj_nz)
    @assert length(adj_I) == length(adj_J) == adj_nnz

    adj_row_sums = adapt(adapter, ci_state.adj_row_sums)
    pos = adapt(adapter, ci_state.pos)
    pos_I = adapt(adapter, ci_state.pos_I)
    pos_J = adapt(adapter, ci_state.pos_J)
    Δpos_IJ = adapt(adapter, ci_state.Δpos_IJ)
    charges = adapt(adapter, ci_state.charges)
    charges_J = adapt(adapter, ci_state.charges_J)
    Tmp_D3_nz = adapt(adapter, ci_state.Tmp_D3_nz)
    S = (x = adapt(adapter, ci_state.S.x), y = adapt(adapter, ci_state.S.y), z = adapt(adapter, ci_state.S.z))
    ϵ_r = adapt(adapter, ci_state.ϵ_r)
    onesT = adapt(adapter, ci_state.onesT)
    contr_thresh = adapt(adapter, ci_state.contr_thresh)
    contr_thresh_I = adapt(adapter, ci_state.contr_thresh_I)

    @inbounds pos_vI = view(pos, adj_I)
    @inbounds pos_vJ = view(pos, adj_J)
    @inbounds charges_vJ = view(charges, adj_J)
    @inbounds contr_thresh_vI = view(contr_thresh, adj_I)

    return ChargeInteractionState(
        M_adj, adj_row_sums, pos, pos_I, pos_J, pos_vI, pos_vJ, Δpos_IJ,
        charges, charges_J, charges_vJ,
        Tmp_D3_nz, S, ϵ_r, onesT, contr_thresh, contr_thresh_I, contr_thresh_vI
    )
end

function resize_charge_interaction_state!!(
    ci_state::ChargeInteractionState{T},
    new_pos::StructVector{<:CartesianPoint{T}}, new_charges::AbstractVector{<:T}, new_ϵ_r::T,
    new_M_adj::AbstractMatrix{T}
) where T
    (;
        adj_row_sums, pos_I, pos_J, Δpos_IJ,
        charges_J,
        Tmp_D3_nz, S, onesT, contr_thresh, contr_thresh_I, contr_thresh_vI
    ) = ci_state

    M_adj = new_M_adj
    pos = new_pos
    charges = new_charges
    ϵ_r = new_ϵ_r

    n_charges = length(pos)

    adj_I, adj_J = findnz(M_adj) # ToDo: Allocation-free change of adj_I and adj_J
    adj_nz = nonzeros(M_adj)
    adj_nnz = length(adj_nz)
    @assert length(adj_I) == length(adj_J) == adj_nnz

    resize!(onesT, n_charges)
    ⋮(onesT) .= one(T)
    resize!(adj_row_sums, n_charges)
    mul!(adj_row_sums, M_adj, onesT)

    resize!(pos_I, adj_nnz)
    @inbounds pos_vI = view(pos, adj_I)
    resize!(pos_J, adj_nnz)
    @inbounds pos_vJ = view(pos, adj_J)
    resize!(Δpos_IJ, adj_nnz)

    resize!(charges_J, adj_nnz)
    @inbounds charges_vJ = view(charges, adj_J)

    resize!(Tmp_D3_nz, adj_nnz)
    S = (
        x = _linked_resize!!(S.x, M_adj),
        y = _linked_resize!!(S.y, M_adj),
        z = _linked_resize!!(S.z, M_adj)
    )

    resize!(contr_thresh, n_charges)
    resize!(contr_thresh_I, adj_nnz)
    @inbounds contr_thresh_vI = view(contr_thresh, adj_I)

    return ChargeInteractionState(
        M_adj, adj_row_sums, pos, pos_I, pos_J, pos_vI, pos_vJ, Δpos_IJ,
        charges, charges_J, charges_vJ,
        Tmp_D3_nz, S, ϵ_r, onesT, contr_thresh, contr_thresh_I, contr_thresh_vI
    )::typeof(ci_state)
end


# ToDo: Make this faster
function dummy_ci_state(pos::StructVector{<:CartesianPoint{T}}, charges::AbstractVector{T}) where T
    # All-zero adjecency matrix, no charge interaction:
    x = pos.x
    M_adj = sparse(Int[], Int[], similar(x, 0), length(x), length(x))
    ϵ_r = one(T)
    return ChargeInteractionState(pos, charges, ϵ_r, M_adj)
end

function full_ci_state(
    pos::StructVector{<:CartesianPoint{T}}, charges::AbstractVector{T}, ϵ_r::T,
) where T
    # Not efficient:
    # M_adj = build_cutoff_matrix(pos, electric_field)

    n_charges = length(charges)
    # ToDo: For charges belonging to multiple events, will need to create a
    # block-diagonal structure:
    M_adj = sparse(fill!(similar(charges, (n_charges, n_charges)), one(T)) - I)

    return ChargeInteractionState(pos, charges, ϵ_r, M_adj)
end


# If field vectors are given, trim adjecency.
# Does not always update charges_J.
function update_charge_interaction!!(
    ci_state::ChargeInteractionState{T},
    new_pos::StructVector{<:CartesianPoint{T}}, new_charges::AbstractVector{<:Real}
) where T
    # ToDo: ignore charges that have been collected (not trapped though!)
    # Will have to remove them from adjecency matrix.

    (;
        M_adj, adj_row_sums, pos, pos_I, pos_J, pos_vI, pos_vJ, Δpos_IJ,
        charges, charges_J, charges_vJ,
        Tmp_D3_nz, S, ϵ_r, onesT, contr_thresh, contr_thresh_I, contr_thresh_vI
    ) = ci_state

    isempty(Tmp_D3_nz) && return ci_state # No interaction, nothing to do.

    parallel_copyto!(pos, new_pos)
    parallel_copyto!(charges, new_charges)

    # Charges may have moved, update pos_I and pos_J from views into charge positions:
    parallel_copyto!(pos_I, pos_vI)
    parallel_copyto!(pos_J, pos_vJ)

    # Update distances between charges:
    adapted_bcast!(-, ⋮, Δpos_IJ.x, pos_I.x, pos_J.x)
    adapted_bcast!(-, ⋮, Δpos_IJ.y, pos_I.y, pos_J.y)
    adapted_bcast!(-, ⋮, Δpos_IJ.z, pos_I.z, pos_J.z)

    inv_4π_ϵ0_ϵ_r = T(inv(4π * ϵ0 * ϵ_r))

    # Update interaction:

    adapted_bcast!(calc_tmp_d3, ⋰, Tmp_D3_nz, Δpos_IJ.x, Δpos_IJ.y, Δpos_IJ.z, inv_4π_ϵ0_ϵ_r)

    adapted_bcast!(*, ⋮, nonzeros(S.x), Tmp_D3_nz, Δpos_IJ.x)
    adapted_bcast!(*, ⋮, nonzeros(S.y), Tmp_D3_nz, Δpos_IJ.y)
    adapted_bcast!(*, ⋮, nonzeros(S.z), Tmp_D3_nz, Δpos_IJ.z)

    return ci_state
end

function calc_tmp_d3(Δx, Δy, Δz, inv_4π_ϵ0_ϵ_r::T) where T
    ϵ_3_2 = T(0.00001)^3  # Set 10µm as the minimum distance
    inv_distance_3 = inv(max((Δx^2 + Δy^2 + Δz^2)^(T(3/2)), ϵ_3_2))
    elementary_charge * inv_4π_ϵ0_ϵ_r * inv_distance_3
end


# Need to call update_charge_interaction!! before and after trim_charge_interaction!!.
function trim_charge_interaction!!(ci_state::ChargeInteractionState{T}, done::AbstractVector{Bool}, field_threshold::T) where T
    (;
        M_adj, adj_row_sums, pos, pos_I, pos_J, pos_vI, pos_vJ, Δpos_IJ,
        charges, charges_J, charges_vJ,
        Tmp_D3_nz, S, ϵ_r, onesT, contr_thresh, contr_thresh_I, contr_thresh_vI
    ) = ci_state

    isempty(nonzeros(M_adj)) && return ci_state # No interaction, nothing to do.

    # Compute thresholds for interaction j -> i:
    # Note: can't use `⋮` for packed Bool vector `done`.
    contr_thresh .= field_threshold ./ adj_row_sums .* _done_field_threshold.(T.(done))
    parallel_copyto!(contr_thresh_I, contr_thresh_vI)
    parallel_copyto!(charges_J, charges_vJ)

    # Trim adjecency matrix:
    nz_S_x, nz_S_y, nz_S_z = ⋮(nonzeros(S.x)), ⋮(nonzeros(S.y)), ⋮(nonzeros(S.z))
    adj_nz = ⋮(nonzeros(M_adj))
    adj_nz .*= (nz_S_x.^2 .+ nz_S_y.^2 .+ nz_S_z.^2) .* ⋮(charges_J).^2 .>= ⋮(contr_thresh_I).^2

    # ToDo: Also use `done` to trim self-interaction sources charges (J indices).

    dropzeros!(M_adj) # ToDo: dropzeros! is not available for CUDA arrays, find alternative.

    return resize_charge_interaction_state!!(ci_state, pos, charges, ϵ_r, M_adj)::typeof(ci_state)
    nothing
end

_done_field_threshold(float_done::Real) = ifelse(iszero(float_done::Real), oftype(float_done::Real, 1), oftype(float_done::Real, Inf))

function apply_charge_interaction!(field_vectors::StructVector{<:CartesianVector{T}}, ci_state::ChargeInteractionState{T}) where T
    (;S, charges) = ci_state

    isempty(nonzeros(S.x)) && return field_vectors # No interaction, nothing to do.

    # In-place muladd, field_vectors.x = S.x * charges * 1 + field_vectors.x * 1, etc:
    mul!(field_vectors.x, S.x, charges, one(T), one(T))
    mul!(field_vectors.y, S.y, charges, one(T), one(T))
    mul!(field_vectors.z, S.z, charges, one(T), one(T))

    return field_vectors
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
    srp_trim_threshold::T = T(20000.0), # 20kV/m
    srp_trim_every::Integer = 100, # trim every 100 steps by default
)::Int where {T <: SSDFloat, S, CC <: ChargeCarrier}
    n_hits::Int, max_nsteps::Int = size(drift_path)
    drift_path[:,1] = startpos
    timestamps[1] = zero(T)
    ϵ_r::T = T(det.semiconductor.material.ϵ_r)

    diffusion_length::T = diffusion ? diff_get_diffusion_length(det, Δt, CC) : zero(T)
    diffusion = diffusion && !iszero(diffusion_length)

    last_real_step_index::Int = 1
    current_pos = StructVector(deepcopy(startpos))
    field_vectors::StructVector{CartesianVector{T}} = StructVector{CartesianVector{T}}(undef, n_hits)
    step_vectors::Vector{CartesianVector{T}} = Vector{CartesianVector{T}}(undef, n_hits)
    done::Vector{Bool} = broadcast(pt -> !_is_next_point_in_det(pt, det, point_types), startpos)
    normal::Vector{Bool} = deepcopy(done)

    if self_repulsion && !all(done)
        ci_state = full_ci_state(current_pos, charges, ϵ_r)
        ci_state = update_charge_interaction!!(ci_state, current_pos, charges)
        ci_state = trim_charge_interaction!!(ci_state, done, srp_trim_threshold)
    else
        # Dummy ChargeInteractionState
        ci_state = dummy_ci_state(current_pos, charges)
    end

    nsteps::Int = 0
    @inbounds for istep in 2:max_nsteps
        nsteps += 1
        last_real_step_index += 1
        fill!(field_vectors, zero(eltype(field_vectors)))
        _add_fieldvector_drift!(field_vectors, current_pos, done, electric_field, det, S)
        if self_repulsion && !all(done)
            # New:
            ci_state = update_charge_interaction!!(ci_state, current_pos, charges)
            apply_charge_interaction!(field_vectors, ci_state)
            if nsteps % srp_trim_every == 0
                ci_state = trim_charge_interaction!!(ci_state, done, srp_trim_threshold)
            end

            # Old:
            # _add_fieldvector_selfrepulsion!(field_vectors, current_pos, done, charges, ϵ_r)
        end

        _get_drift_steps!(step_vectors, field_vectors, done, Δt, det.semiconductor.charge_drift_model, CC)
        diffusion && _add_fieldvector_diffusion!(step_vectors, done, diffusion_length)
        _modulate_driftvectors!(step_vectors, current_pos, det.virtual_drift_volumes)
        _check_and_update_position!(step_vectors, current_pos, done, normal, drift_path, timestamps, istep, det, grid, point_types, startpos, Δt, verbose)
        if all(done) break end
    end

    return last_real_step_index
end


#!!!!!!!!!!!!!! REMOVE THIS
function _drift_charge_oldsrp!(
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
)::Int where {T <: SSDFloat, S, CC <: ChargeCarrier}
    n_hits::Int, max_nsteps::Int = size(drift_path)
    drift_path[:,1] = startpos
    timestamps[1] = zero(T)
    ϵ_r::T = T(det.semiconductor.material.ϵ_r)

    diffusion_length::T = diffusion ? diff_get_diffusion_length(det, Δt, CC) : zero(T)
    diffusion = diffusion && !iszero(diffusion_length)

    last_real_step_index::Int = 1
    current_pos::Vector{CartesianPoint{T}} = deepcopy(startpos)
    field_vectors::Vector{CartesianVector{T}} = Vector{CartesianVector{T}}(undef, n_hits)
    step_vectors::Vector{CartesianVector{T}} = Vector{CartesianVector{T}}(undef, n_hits)
    done::Vector{Bool} = broadcast(pt -> !_is_next_point_in_det(pt, det, point_types), startpos)
    normal::Vector{Bool} = deepcopy(done)

    @inbounds for istep in 2:max_nsteps
        last_real_step_index += 1
        fill!(field_vectors, zero(eltype(field_vectors)))
        _add_fieldvector_drift!(field_vectors, current_pos, done, electric_field, det, S)
        self_repulsion && _add_fieldvector_selfrepulsion!(field_vectors, current_pos, done, charges, ϵ_r)
        _get_drift_steps!(step_vectors, field_vectors, done, Δt, det.semiconductor.charge_drift_model, CC)
        diffusion && _add_fieldvector_diffusion!(step_vectors, done, diffusion_length)
        _modulate_driftvectors!(step_vectors, current_pos, det.virtual_drift_volumes)
        _check_and_update_position!(step_vectors, current_pos, done, normal, drift_path, timestamps, istep, det, grid, point_types, startpos, Δt, verbose)
        if all(done) break end
    end

    return last_real_step_index
end

function _add_fieldvector_drift!(step_vectors::AbstractVector{CartesianVector{T}}, current_pos::AbstractVector{CartesianPoint{T}}, done::AbstractVector{Bool}, electric_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, det::SolidStateDetector{T}, ::Type{S})::Nothing where {T, S}
    #for n in eachindex(step_vectors)
    #   if !done[n]
    #       step_vectors[n] += get_velocity_vector(electric_field, _convert_point(current_pos[n], S))
    #       done[n] = (step_vectors[n] == CartesianVector{T}(0,0,0))
    #   end
    #end
    step_vectors .= step_vectors .+ get_velocity_vector_withconv.(Ref(electric_field), current_pos, S, done)
    done .= iszero.(step_vectors)
    nothing
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
