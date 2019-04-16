# function drift_charge!(
#                             drift_path::Vector{CartesianPoint{T}},
#                             det::SolidStateDetector{T, :Cylindrical},
#                             startpos::CartesianPoint{T},
#                             delta_t::T,
#                             velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}
#                         )::Nothing where {T <: SSDFloat}
#     # drifttime::T = 0.0
#     done::Bool = false
#     stepvector = @SVector zeros(T, 3)
#     pos_cyl = CylindricalPoint{T}(0, 0, 0)
#     crossing_pos_cyl = CylindricalPoint{T}(0, 0, 0)
#     crossing_pos_type::Symbol = :Idle
#     boundary_index::Int=0
#     stepvector_ref = SVector{3,T}(0.0,0.0,0.0)
#     drift_path[1] = startpos
#     @show startpos
#     for istep in eachindex(drift_path)[2:end]
#         if done == false
#             pos_cyl = CylindricalPoint(CartesianPoint{T}(drift_path[istep - 1])) # update pos_cyl
#             # pos_cyl = geom_round(CylindricalPoint(CartesianPoint{T}(drift_path[istep-1]))) # update pos_cyl
#             if pos_cyl in det
#                 stepvector = get_velocity_vector(velocity_field, pos_cyl) * delta_t
#                 drift_path[istep] = drift_path[istep-1] + stepvector
#                 stepvector == stepvector_ref ? done = true : nothing
#                 # stepvector == stepvector_ref ? drifttime = delta_t * (istep-1) : nothing
#                 # stepvector == stepvector_ref ? println(pos_cyl) : nothing
#             elseif istep < 3
#                 done = true
#                 # drifttime = delta_t * (istep-1)
#                 throw(ErrorException("start position $pos_cyl is not contained in detector"))
#             else
#                 crossing_pos_cyl, crossing_pos_type, boundary_index = get_crossing_pos(det, drift_path[istep-2], drift_path[istep-1])
#                 # println(crossing_pos_cyl, crossing_pos_type, boundary_index)
#                 if crossing_pos_type == :electrode
#                     done = true
#                     # drifttime = delta_t * (istep-1)
#                     drift_path[istep-1] = CartesianPoint(crossing_pos_cyl)
#                     drift_path[istep] = drift_path[istep-1]
#                 elseif crossing_pos_type == :floating_boundary
#                     # println("outch")
#                     # println(crossing_pos_cyl)
#                     stepvector = get_velocity_vector(velocity_field, crossing_pos_cyl) * delta_t
#                     # stepvector = get_velocity_vector(velocity_field, CylindricalPoint{T}(drift_path[istep-2]...)) * delta_t
#                     # projected_vector = project_to_plane(stepvector, get_fb_plane(crossing_pos_cyl, drift_path[istep-2], drift_path[istep-1]).n⃗)
#                     projected_vector = project_to_plane(stepvector, is_surface_point(det, crossing_pos_cyl)[2] )
#                     # drift_path[istep-1] = CartesianPoint(crossing_pos_cyl)
#                     # drift_path[istep] = drift_path[istep-1] + projected_vector

#                     drift_path[istep-1] = drift_path[istep-2]+ projected_vector
#                     increment::T=2.0
#                     while !(in(CylindricalPoint(CartesianPoint{T}(drift_path[istep-1]...)), det)) && increment < 500
#                         # println(crossing_pos_cyl)
#                         drift_path[istep-1] = drift_path[istep-2]+ 1.0/increment * projected_vector
#                         # println(drift_path[istep-1])
#                         increment+=1
#                     end
#                     increment == 500 ? @warn("Handling of charge at floating boundary did not work as intended. Start Position (Cart): $startpos") : nothing
#                     drift_path[istep] = drift_path[istep-1] + projected_vector
#                 else
#                     @show crossing_pos_cyl,crossing_pos_type, boundary_index
#                     @show typeof(crossing_pos_cyl)
#                     @show typeof(pos_cyl)
#                     @warn ("Internal error for charge stating at $startpos")
#                     drift_path[istep] = drift_path[istep-1]
#                     # drifttime = delta_t * (istep-1)
#                     done = true
#                 end
#             end
#         elseif done == true
#             drift_path[istep] = drift_path[istep-1]
#         end
#     end
#     # drifttime
#     return nothing
# end

"""
    Before calling this function one should check that `startpos` is inside `det`: `in(startpos, det`

"""
function drift_charge!(
                            drift_path::Vector{CartesianPoint{T}},
                            det::SolidStateDetector{T, :Cylindrical},
                            grid::Grid{T, 3, :Cylindrical},
                            startpos::CartesianPoint{T},
                            Δt::T,
                            velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}
                        )::Nothing where {T <: SSDFloat}
    drifttime::T = 0
    done::Bool = false
    drift_path[1] = startpos
    null_step::CartesianVector{T} = CartesianVector{T}(0, 0, 0)
    for istep in eachindex(drift_path)[2:end] #end] 
        if done == false
            current_pos::CartesianPoint{T} = drift_path[istep - 1]
            stepvector::CartesianVector{T} = get_velocity_vector(velocity_field, CylindricalPoint(current_pos)) * Δt
            if stepvector == null_step done = true end 
            next_pos::CartesianPoint{T} = current_pos + stepvector
            if CylindricalPoint(next_pos) in det
                drift_path[istep] = next_pos
            else
                crossing_pos::CartesianPoint{T}, cd_point_type::UInt8, boundary_index::Int, surface_normal::CartesianVector{T} = get_crossing_pos(det, grid, current_pos, next_pos)
                if cd_point_type == CD_ELECTRODE
                    done = true
                    drift_path[istep] = crossing_pos
                elseif cd_point_type == CD_FLOATING_BOUNDARY
                    projected_vector = CartesianVector{T}(project_to_plane(stepvector, surface_normal))
                    next_pos = current_pos + projected_vector
                    # ToDo: We actually need a time array as well to do this properly...
                    small_projected_vector = projected_vector * T(0.001)
                    i::Int = 0
                    while !(CylindricalPoint(next_pos) in det) && i < 1000
                        next_pos -= small_projected_vector
                        i += 1
                    end
                    if i == 1000 @warn("Handling of charge at floating boundary did not work as intended. Start Position (Cart): $startpos") end
                    drift_path[istep] = next_pos
                else 
                    @warn ("Internal error for charge stating at $startpos")
                    done = true
                end
            end
        elseif done == true
            drift_path[istep] = drift_path[istep-1]
        end
    end
    return nothing
end

# function get_crossing_pos(  detector::SolidStateDetector{T, :Cylindrical}, point_in::CartesianPoint{T}, point_out::CartesianPoint{T}; 
#                             dig::Int=6, max_n_iter::Int=2000) where {T <: SSDFloat}
#     T==Float64 ? dig = 13 : nothing
#     current_pos_in = MVector{3,T}(point_in...)
#     current_pos_out = MVector{3,T}(point_out...)
#     current_pos_xyz::SVector{3,T} = 0.5 * (current_pos_in .+ current_pos_out)
#     current_pos_cyl::CylindricalPoint{T} = CylindricalPoint(CartesianPoint{T}(current_pos_xyz))
#     i_limit::Int = 0
#     while (point_type(detector, current_pos_cyl)[2] == -1) && i_limit < max_n_iter
#         if contains(detector, current_pos_cyl)
#             current_pos_in[1] = current_pos_xyz[1]
#             current_pos_in[2] = current_pos_xyz[2]
#             current_pos_in[3] = current_pos_xyz[3]
#         else
#             current_pos_out[1] = current_pos_xyz[1]
#             current_pos_out[2] = current_pos_xyz[2]
#             current_pos_out[3] = current_pos_xyz[3]
#         end
#         current_pos_xyz = 0.5 * (current_pos_out .+ current_pos_in)
#         current_pos_cyl = (CylindricalPoint(CartesianPoint{T}(current_pos_xyz)))

#         i_limit += 1
#     end
#     dig=6
#     crossing_pos_type, boundary_index = point_type(detector,current_pos_cyl)
#     if i_limit==max_n_iter
#         @warn "kek"
#         current_pos_cyl, crossing_pos_type, boundary_index = get_crossing_pos_w_rounding(detector, point_in, point_out)
#     end
#     # geom_round(current_pos_cyl), crossing_pos_type, boundary_index
#     return current_pos_cyl, crossing_pos_type, boundary_index
# end

# Point types for charge drift: Defined in DetectorGeometries/DetectorGeometries.jl
# const CD_ELECTRODE = 0x00
# const CD_OUTSIDE = 0x01
# const CD_BULK = 0x02
# const CD_FLOATING_BOUNDARY = 0x04 # not 0x03, so that one could use bit operations here...

function get_crossing_pos(  detector::SolidStateDetector{T, :Cylindrical}, grid::Grid{T, 3}, point_in::CartesianPoint{T}, point_out::CartesianPoint{T}; 
                            max_n_iter::Int = 2000)::Tuple{CartesianPoint{T}, UInt8, Int, CartesianVector{T}} where {T <: SSDFloat}
    point_mid::CartesianPoint{T} = T(0.5) * (point_in + point_out)
    cd_point_type::UInt8, contact_idx::Int, surface_normal::CartesianVector{T} = point_type(detector, grid, CylindricalPoint(point_mid))
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
        cd_point_type, contact_idx, surface_normal = point_type(detector, grid, CylindricalPoint(point_mid))
    end
    return point_mid, cd_point_type, contact_idx, surface_normal
end



function drift_and_acc_charge!(
    drift_path::Vector{CartesianPoint{T}},
    charge_signal::AbstractVector{<:RealQuantity},
    detector::SolidStateDetector,
    velocity_field::Interpolations.Extrapolation{<:SVector{3},3},
    weighting_potential::Interpolations.Extrapolation{<:Real,3},
    startpos::CartesianPoint{T},
    charge::RealQuantity,
    delta_t::RealQuantity,
) where {T <: SSDFloat}
    float_delta_t = float(delta_t)
    time_unit = unit(float_delta_t)
    unitless_delta_t = to_internal_units(u"s", float_delta_t)

    unitless_hitpos = to_internal_units.(u"m", startpos)

    t_drift = drift_charge!(
        drift_path, detector, unitless_hitpos,
        unitless_delta_t, velocity_field
    )

    accumulate_charge!(charge_signal, drift_path, weighting_potential, charge)

    from_internal_units(time_unit, u"s", t_drift)
end



struct ChargeDriftEvent{T<:Real}
    drift_paths::AbstractVector{<:DriftPath}
    signal::AbstractVector{<:AbstractVector{T}}

    ChargeDriftEvent{T}(drift_paths::AbstractVector{<:DriftPath},signal::AbstractVector{<:AbstractVector{T}}) where {T<:Real} = new{T}(drift_paths,signal)
end

function ChargeDriftEvent(drift_paths::AbstractVector{<:DriftPath},signal::AbstractVector{<:AbstractVector{T}}) where T<:Real
    return ChargeDriftEvent{T}(drift_paths, signal)
end

function pulse_from_drift_paths(drift_paths, energy_depositions::AbstractVector{T}, Wpot_interp::Interpolations.Extrapolation{T,3}, return_contributions::Bool = false) where T <:Real
    if return_contributions == true
        charge_signal_e = zeros(T,size(drift_paths[1].e_path,1))
        charge_signal_h = zeros(T,size(drift_paths[1].e_path,1))
        for i in eachindex(drift_paths)
            accumulate_charge!(charge_signal_e, drift_paths[i].e_path, Wpot_interp, energy_depositions[i])
            accumulate_charge!(charge_signal_h, drift_paths[i].h_path, Wpot_interp, energy_depositions[i])
        end
        # charge_signal = charge_signal_h .-= charge_signal_e
        return charge_signal_e, charge_signal_h
    else
        charge_signal = zeros(T,size(drift_paths[1].e_path,1))
        for i in eachindex(drift_paths)
            accumulate_charge!(charge_signal, drift_paths[i].e_path, Wpot_interp, -1*energy_depositions[i])
            accumulate_charge!(charge_signal, drift_paths[i].h_path, Wpot_interp, energy_depositions[i])
        end
        return charge_signal
    end
end


# function drift_charges( detector::SolidStateDetector{T}, starting_positions::Vector{CylindricalPoint{T}}, 
#                         velocity_field_e::Interpolations.Extrapolation{SVector{3, Float64}, 3},
#                         velocity_field_h::Interpolations.Extrapolation{SVector{3, Float64}, 3}; 
#                         Δt::T = T(1f-9), n_steps::Int = 2000) where {T <: SSDFloat}

#     drift_path_e = [@SVector zeros(T,3) for i in 1:n_steps]
#     drift_path_h = [@SVector zeros(T,3) for i in 1:n_steps]

#     for j in eachindex(starting_positions)
#         drift_charge!(drift_path_e, detector, starting_positions[j], delta_t, velocity_field_e)
#         drift_charge!(drift_path_h, detector, starting_positions[j], delta_t, velocity_field_h)
#         drift_paths[j] = DriftPath(deepcopy(drift_path_e), deepcopy(drift_path_h))
#     end
#     drift_paths
# end


"""
    ToDo: Change name to something like get_signals(W_pots, drift_path)...
"""
function drift_charges( detector::SolidStateDetector{T}, 
                        grid::Grid{T, 3},
                        starting_positions::Vector{CylindricalPoint{T}}, 
                        velocity_field_e::Interpolations.Extrapolation{SVector{3,Float64},3}, velocity_field_h::Interpolations.Extrapolation{SVector{3,Float64},3},
                        energy_depositions::AbstractVector{T}, weighting_potentials::AbstractVector{<:Interpolations.Extrapolation{T,3}}; 
                        delta_t::T = T(1f-9), n_steps::Int = 2000) where {T <: SSDFloat}


    drift_paths = drift_charges(detector, grid, starting_positions, velocity_field_e, velocity_field_h, delta_t=delta_t, n_steps=n_steps)
    signal=Vector{AbstractVector{T}}(undef,size(weighting_potentials,1))
    for i in eachindex(weighting_potentials)
        signal[i]=pulse_from_drift_paths(drift_paths, energy_depositions, weighting_potentials[i])
    end
    return ChargeDriftEvent(drift_paths,signal)
end

function get_crossing_pos_w_rounding(detector::SolidStateDetector, grid::Grid{T, 3}, point_in::CartesianPoint{T}, point_out::CartesianPoint{T}; max_n_iter::Int=2000) where T <:Real
    current_pos_in = MVector{3,T}(point_in...)
    current_pos_out = MVector{3,T}(point_out...)
    current_pos_xyz::SVector{3,T} = 0.5 * (current_pos_in .+ current_pos_out)
    current_pos_cyl::CylindricalPoint{T} = geom_round(CylindricalPoint(CartesianPoint{T}(current_pos_xyz)))
    i_limit::Int=0
    println("intersection w rounding...")
    while (point_type(detector, grid, current_pos_cyl)[2]==-1) && i_limit < max_n_iter
        if contains(detector,current_pos_cyl)
            current_pos_in[1] = current_pos_xyz[1]
            current_pos_in[2] = current_pos_xyz[2]
            current_pos_in[3] = current_pos_xyz[3]
        else
            current_pos_out[1] = current_pos_xyz[1]
            current_pos_out[2] = current_pos_xyz[2]
            current_pos_out[3] = current_pos_xyz[3]
        end
        current_pos_xyz = 0.5 * (current_pos_out .+ current_pos_in)
        current_pos_cyl=geom_round(CylindricalPoint(CartesianPoint{T}(current_pos_xyz)))
        i_limit +=1
    end
    if i_limit == max_n_iter
        @warn("charge cannot be collected, something is wrong...")
    end
    crossing_pos_type, boundary_index = point_type(detector, grid, current_pos_cyl)
    println(current_pos_xyz)
    current_pos_cyl, crossing_pos_type, boundary_index
end

function project_to_plane(v⃗::AbstractArray, n⃗::AbstractArray) #Vector to be projected, #normal vector of plane
    # solve (v⃗+λ*n⃗) ⋅ n⃗ = 0
    # n⃗ = n⃗ ./ norm(n⃗)
    λ=-1*dot(v⃗,n⃗)/dot(n⃗,n⃗)
    SVector{3,eltype(v⃗)}(v⃗[1]+λ*n⃗[1],v⃗[2]+λ*n⃗[2],v⃗[3]+λ*n⃗[3])
end

function lineplanecollision(planenorm::Vector, planepnt::Vector, raydir::Vector, raypnt::Vector)
    ndotu = dot(planenorm, raydir)
    if ndotu ≈ 0 error("no intersection or line is within plane") end

    w  = raypnt - planepnt
    si = -dot(planenorm, w) / ndotu
    ψ  = w .+ si .* raydir .+ planepnt
    return ψ
end

# function get_cartesian_paths(a::Array{CylindricalPoint,1})
#     xpath = Array{eltype(a[1].r),1}(undef,size(a))
#     ypath = Array{eltype(a[1].r),1}(undef,size(a))
#     zpath = Array{eltype(a[1].r),1}(undef,size(a))
#     xyz = map(x->CartesianPoint(x),a)
#     for i in eachindex(xyz)
#         xpath[i]=xyz[i].x
#         ypath[i]=xyz[i].y
#         zpath[i]=xyz[i].z
#     end
#     xpath,ypath,zpath
# end
# function get_cartesian_paths(a::Array{SArray{Tuple{3},T,1,3},1};scaling::Real = 1) where T <:Real
#     xpath = Array{T,1}(undef,size(a))
#     ypath = Array{T,1}(undef,size(a))
#     zpath = Array{T,1}(undef,size(a))
#     for i in eachindex(a)
#         xpath[i]=a[i][1]*scaling
#         ypath[i]=a[i][2]*scaling
#         zpath[i]=a[i][3]*scaling
#     end
#     xpath,ypath,zpath
# end
function get_rφz_vector_from_xyz(vector::AbstractArray)::AbstractArray
  @fastmath begin
    x = vector[1]
    y = vector[2]
    z = vector[3]
    r = sqrt(x^2 + y^2)
    # φ = atan(y , x)
    φ = atan(y / x)
    x <  0 && y >= 0 ? φ += π  : nothing
    x <= 0 && y <  0 ? φ += π  : nothing
    x >  0 && y <  0 ? φ += 2π : nothing
  end
  return [r, φ, z]
end

function get_rφz_from_xyz(v::AbstractArray)::AbstractArray
    r = sqrt(v[1]^2+v[2]^2)
    φ = atan(v[2],v[1])
    while φ<0
        φ+=2π
    end
    while φ>2π
        φ-=2π
    end
    z=v[3]
    [r,φ,z]
end

function convert_xyz_vectorfield_in_rφz_to_rφz_vectorfield_in_rφz(field,grid)
    outputfield = Array{Array{Float32,1}}(undef,size(field,1),size(field,2),size(field,3))
    for (iz,z) in enumerate(grid.z)
        for (iφ,φ) in enumerate(grid.φ)
            for (ir,r) in enumerate(grid.r)
                outputfield[ir,iφ,iz]=get_rφz_from_xyz(field[ir,iφ,iz])
            end
        end
    end
    outputfield
end

function rotate_vector_z_axis(vector,α)
    Rα = Array{AbstractFloat}(undef,2,2)
    Rα[1,1]=cos(α)
    Rα[1,2]=-1*sin(α)
    Rα[2,1]=sin(α)
    Rα[2,2]=cos(α)
    result = Rα*vector[1:2]
    push!(result,vector[3])
    result
end

function add_cylindric_vectors(v1,v2)::AbstractArray
    r = sqrt(v1[1]^2 +v2[1]^2 + 2*v1[1]*v2[1]*cos(v2[2]-v1[2]))
    φ = v1[2]+atan(v2[1]*sin(v2[2]-v1[2]),v1[1]+v2[1]*cos(v2[2]-v1[2]))
    z = v1[3]+v2[3]
    [r,φ,z]
end


function get_velocity_vector(interpolation_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, point::CartesianPoint{T})::CartesianVector{T} where {T <: SSDFloat}  
    return CartesianVector{T}(interpolation_field(point.x, point.y, point.z))
end

@inline get_velocity_vector(interpolated_vectorfield, point::CylindricalPoint{T}) where {T <: SSDFloat} =
    CartesianVector{T}(interpolated_vectorfield(point.r, point.φ, point.z))


struct Trajectory
    startposition::Vector
    n_charges::Int
    indv_charge::Int
    n_steps::Int
    timestep::AbstractFloat
    rφz_path::Array{CylindricalPoint,1}
    x_path::Vector
    y_path::Vector
    z_path::Vector
    endposition::Vector
    drifttime::AbstractFloat

    function Trajectory(startposition, n_charges, indv_charge, n_steps, timestep, rφz_path, x_path, y_path, z_path,endposition, drifttime)
        return new(startposition, n_charges, indv_charge, n_steps, timestep, rφz_path, x_path,y_path, z_path, endposition, drifttime)
    end
end

function println(io::IO, trj::Trajectory)
    println("____Trajectory_____")
    # print("Drift of $(trj.n_charges) ")
    trj.indv_charge == 1 ? println("h+") : println("e-")
    println("Startposition: x:$(round(trj.startposition[1],digits=5)) y:$(round(trj.startposition[2],digits=5)) z:$(round(trj.startposition[3],digits=5))")
    println("Endposition: x:$(round(trj.endposition[1],digits=5)) y:$(round(trj.endposition[2],digits=5)) z:$(round(trj.endposition[3],digits=5))")
    println("Drifttime: $(trj.drifttime) ns")
end
function show(io::IO, trj::Trajectory) println(trj) end
function print(io::IO, trj::Trajectory) println(trj) end
function display(io::IO, trj::Trajectory) println(trj) end
function show(io::IO,::MIME"text/plain", trj::Trajectory)
    show(io, trj)
end

struct Pulse
    electron_contribution::Vector
    hole_contribution::Vector
    signal::Vector
    energy::AbstractFloat
    function Pulse(electron_contribution, hole_contribution, signal,energy)
        return new(electron_contribution, hole_contribution, signal,energy)
    end
end

@inline function _get_wpot_at(Wpot_interp::Interpolations.Extrapolation{T,3}, pos::SVector{3,T}) where T<:Real
    pos_cyl::CylindricalPoint{T} = CylindricalPoint(CartesianPoint{T}(pos)) #CylFromCart(pos)
    # isapprox(Wpot_interp(pos_cyl.r,pos_cyl.φ,pos_cyl.z),0.0,atol = 5e-5) ? T(0.0) : Wpot_interp(pos_cyl.r,pos_cyl.φ,pos_cyl.z)
    Wpot_interp(pos_cyl.r,pos_cyl.φ,pos_cyl.z)
end


function accumulate_charge!(charge_signal::AbstractVector{<:RealQuantity}, drift_path::AbstractVector{<:SVector{3,T}}, Wpot_interp::Interpolations.Extrapolation{T,3}, charge::RealQuantity) where T <:Real
    idxs = eachindex(charge_signal)
    eachindex(drift_path) == idxs || throw(ArgumentError("Vectors charge_signal and drift_path must have same indices"))
    pos_prev = drift_path[first(idxs)]
    V_0 = _get_wpot_at(Wpot_interp, pos_prev)
    V_prev::typeof(V_0) = V_0
    @inbounds for i in idxs
        pos_current = drift_path[i]
        V_current = _get_wpot_at(Wpot_interp, drift_path[i])
        charge_signal[i] = muladd(charge, V_current - V_0, charge_signal[i])
    end
    charge_signal
end



function generate_charge_signals!(
    contact_charge_signals::AbstractVector{<:AbstractVector{<:AbstractVector{<:RealQuantity}}},
    drift_times::AbstractVector{<:RealQuantity},
    detector::SolidStateDetector,
    velocity_field_e::Interpolations.Extrapolation{<:SVector{3},3},
    velocity_field_h::Interpolations.Extrapolation{<:SVector{3},3},
    weighting_potentials::AbstractVector{<:Interpolations.Extrapolation{<:Real,3}},
    hit_pos::Vector{Vector{CartesianPoint{T}}},
    hit_edep::AbstractVector{<:AbstractVector{<:RealQuantity}},
    E_ionisation::RealQuantity,
    delta_t::RealQuantity,
) where {T <: SSDFloat}
    float_delta_t = float(delta_t)
    time_unit = unit(float_delta_t)
    unitless_delta_t = to_internal_units(u"s", float_delta_t)

    T_len = eltype(eltype(velocity_field_e.itp.coefs))
    drift_path = Vector{SVector{3,T_len}}()

    for i_evt in eachindex(hit_pos)
        first_contact_idx = first(eachindex(weighting_potentials))
        resize!(drift_path, size(contact_charge_signals[first_contact_idx][i_evt], 1))

        t_drift_total::eltype(drift_times) = zero(eltype(drift_times))
        timeunit = unit(t_drift_total)

        for i_wpot in eachindex(first_contact_idx)
            fill!(contact_charge_signals[i_wpot][i_evt], 0)
        end

        startpos_vec = hit_pos[i_evt]
        edep_vec = hit_edep[i_evt]
        start_positions = hit_pos[i_evt]
        for i_hit in eachindex(start_positions)
            unitless_startpos = CartesianPoint{Float32}(to_internal_units(u"m", startpos_vec[i_hit]))
            edep = edep_vec[i_hit]
            n_charges = uconvert(NoUnits,edep / E_ionisation)

            t_drift_e = 1u"s" * drift_charge!(
                drift_path, detector, unitless_startpos,
                unitless_delta_t, velocity_field_e
            )
            for i_wpot in eachindex(weighting_potentials)
                signal = contact_charge_signals[i_wpot][i_evt]
                accumulate_charge!(signal, drift_path, weighting_potentials[i_wpot], -n_charges)
            end
            t_drift_total = max(t_drift_total, uconvert(timeunit, t_drift_e))

            t_drift_h = 1u"s" * drift_charge!(
                drift_path, detector, unitless_startpos,
                unitless_delta_t, velocity_field_h
            )
            for i_wpot in eachindex(weighting_potentials)
                signal = contact_charge_signals[i_wpot][i_evt]
                accumulate_charge!(signal, drift_path, weighting_potentials[i_wpot], +n_charges)
            end
            t_drift_total = max(t_drift_total, uconvert(timeunit, t_drift_h))
        end

        drift_times[i_evt] = t_drift_total
    end
    nothing
end


function generate_charge_signals(
    detector::SolidStateDetector,
    velocity_field_e::Interpolations.Extrapolation{<:SVector{3},3},
    velocity_field_h::Interpolations.Extrapolation{<:SVector{3},3},
    weighting_potentials::AbstractVector{<:Interpolations.Extrapolation{<:Real,3}},
    events::DetectorHitEvents,
    n_steps::Integer,
    delta_t::RealQuantity,
)
    T_charge = Float32
    T_time = Float32

    E_ionisation = detector.material_detector.E_ionisation

    hit_pos = events.pos
    hit_edep = events.edep

    drift_times = fill(zero(T_time) * u"ns", size(hit_pos, 1))
    contact_charge_signals = [nestedview(Array{T_charge}(undef, n_steps, size(hit_pos, 1))) for i in eachindex(weighting_potentials)]

    generate_charge_signals!(
        contact_charge_signals, drift_times,
        detector, velocity_field_e, velocity_field_h, weighting_potentials,
        hit_pos, hit_edep,
        E_ionisation, delta_t
    )

    contact_charge_signals, drift_times
end


function Pulse(Trajectory_e::Trajectory, Trajectory_h::Trajectory, wp_int, energy::T = 1.0) where {T <: SSDFloat}#wp_int = interpolation field for the weighting potential
    n_steps = Trajectory_e.n_steps
    electron_contribution::Array{T,1} = zeros(T, n_steps)
    hole_contribution::Array{T,1} = zeros(T, n_steps)
    signal::Array{T,1} = zeros(T,n_steps)
    for i_step in eachindex(electron_contribution)
        curr_loc_xyz_e::SVector{3,T} = SVector(Trajectory_e.x_path[i_step],Trajectory_e.y_path[i_step],Trajectory_e.z_path[i_step])
        # curr_loc_rφz_e=get_rφz_from_xyz(curr_loc_xyz_e)
        curr_loc_rφz_e::CylindricalPoint{T} = CylindricalPoint(CartesianPoint{T}(curr_loc_xyz_e)) #CylFromCart(curr_loc_xyz_e)
        curr_loc_xyz_h::SVector{3,T} = SVector(Trajectory_h.x_path[i_step],Trajectory_h.y_path[i_step],Trajectory_h.z_path[i_step])
        # curr_loc_rφz_h=get_rφz_from_xyz(curr_loc_xyz_h)
        curr_loc_rφz_h::CylindricalPoint{T} = CylindricalPoint(CartesianPoint{T}(curr_loc_xyz_h)) #CylFromCart(curr_loc_xyz_h)
        # electron_contribution[i_step] = -1 * energy * wp_int(curr_loc_rφz_e...)
        electron_contribution[i_step] = -1 * energy * wp_int(curr_loc_rφz_e.r, curr_loc_rφz_e.φ, curr_loc_rφz_e.z)
        hole_contribution[i_step] = energy * wp_int(curr_loc_rφz_h.r, curr_loc_rφz_h.φ, curr_loc_rφz_h.z)
        signal[i_step] = (hole_contribution[i_step] + electron_contribution[i_step])
    end
    Pulse(SVector{n_steps}(electron_contribution), SVector{n_steps}(hole_contribution), SVector{n_steps}(signal), energy)
end

function println(io::IO,P::Pulse)
    println("test")
end


struct Energy_Deposition
    location::SVector{3}
    energy::AbstractFloat
    n_eh_pairs::Int
    function Energy_Deposition(location::AbstractArray,energy::AbstractFloat,n_eh_pairs::Int)
        return new(location,energy,n_eh_pairs)
    end
    function Energy_Deposition(location::CylindricalPoint{T},energy::AbstractFloat,n_eh_pairs::Int) where {T <: SSDFloat}
        return new(CartesianPoint(location),energy,n_eh_pairs)
    end
end

function println(io::IO,edep::Energy_Deposition)
    println("___Energy_Deposition___")
    println("Position: x:$(round(edep.location[1], digits=6)) y:$(round(edep.location[2], digits=6)) z:$(round(edep.location[3], digits=6))")
    println("Energy: $(edep.energy) keV")
end
function show(io::IO, edep::Energy_Deposition) println(edep) end
function print(io::IO, edep::Energy_Deposition) println(edep) end
function display(io::IO, edep::Energy_Deposition) println(edep) end
function show(io::IO,::MIME"text/plain", edep::Energy_Deposition)
    show(io, edep)
end

mutable struct Event
    detector::SolidStateDetector
    n_sites::Int
    energy_depositions::Array{Energy_Deposition,1}
    total_energy::AbstractFloat
    velocity_field_e
    velocity_field_h
    trajectories_e::Array{Trajectory,1}
    trajectories_h::Array{Trajectory,1}
    weighting_potentials::AbstractArray
    pulses::Array{Pulse,1}
    function Event()
        return new()
    end
end

function println(io::IO, E::Event)
    println("_________Event___________")
    println("Energy:\t$(round(E.total_energy,digits=5)) keV")
    println("RealQuantity of Hits:\t$(E.n_sites)")
end
function show(io::IO, E::Event) println(E) end
function print(io::IO, E::Event) println(E) end
function display(io::IO, E::Event) println(E) end
function show(io::IO,::MIME"text/plain", E::Event)
    show(io, E)
end

function +(p1::Pulse,p2::Pulse)
    return Pulse(p1.electron_contribution .+ p2.electron_contribution, p1.hole_contribution .+ p2.hole_contribution, p1.signal .+ p2.signal, p1.energy + p2.energy)
end

function Event(det::SolidStateDetector, energy_depositions::Array{Energy_Deposition,1}, int_velocity_field_e, int_velocity_field_h, weighting_potentials; n_steps::Int=1000,timestep=1e-9)
    e = Event()
    e.detector=det
    e.trajectories_e =[]
    e.trajectories_h =[]
    e.pulses = []
    e.energy_depositions = energy_depositions
    e.total_energy = sum(map(x->x.energy,energy_depositions))
    e.n_sites = size(energy_depositions,1)
    e.velocity_field_e = int_velocity_field_e
    e.velocity_field_h = int_velocity_field_h
    for it in 1:e.n_sites
        push!(e.trajectories_e,Trajectory(det,energy_depositions[it].location,n_steps,timestep,:e,int_velocity_field_e))
        push!(e.trajectories_h,Trajectory(det,energy_depositions[it].location,n_steps,timestep,:h,int_velocity_field_h))
    end

    e.weighting_potentials = weighting_potentials
    for iseg in eachindex(weighting_potentials)
        trj_pulses = []
        for itr in 1:e.n_sites
            push!(trj_pulses,Pulse(e.trajectories_e[itr],e.trajectories_h[itr],e.weighting_potentials[iseg],e.energy_depositions[itr].energy))
        end
        push!(e.pulses,sum(trj_pulses))
    end
    e
end

function Event(det::SolidStateDetector, starting_positions::AbstractArray, energies, int_velocity_field_e, int_velocity_field_h, weighting_potentials; n_steps=1000, timestep=1e-9)
    energy_depositions::Array{Energy_Deposition,1} = []
    for i in eachindex(starting_positions)
        push!(energy_depositions, Energy_Deposition(starting_positions[i], energies[i], round(Int, energies[i]/3.0e-3)))
    end
    Event(det, energy_depositions, int_velocity_field_e, int_velocity_field_h, weighting_potentials,n_steps=n_steps,timestep=timestep)
end

struct StraightLine
    reference_point
    step_vector
    function StraightLine(reference_point,step_vector)
        return new(reference_point,step_vector)
    end
end
function straight_line_from_points(p1,p2)
    reference_point = p1
    step_vector = (p2-p1)
    step_vector=step_vector ./ norm(step_vector)
    StraightLine(reference_point,step_vector)
end

struct Plane
    reference_point::AbstractArray
    n⃗::AbstractArray
    function Plane(n⃗,reference_point)
        return new(n⃗,reference_point)
    end
end

function get_fb_plane(pCyl::CylindricalPoint{T}, p_in::SVector{3,T}, p_out::SVector{3,T})::Plane where {T <: SSDFloat}
    pCart::CartesianPoint{T} = CartesianPoint(pCyl)
    return Plane(pCart,SVector{3,Float32}(normalize!([[p_out-p_in]...])...))
end

function cramers_rule_2d(M,V)
    nominator = det(M)
    Mx = hcat(V,M[:,2])
    My = hcat(M[:,1],V)
    result= @MVector [det(Mx), det(My)]
    result = result ./ nominator
end

function cramers_rule_3d(M,V)
    nominator = det(M)
    Mx = hcat(V',M[:,2:3])
    My = hcat(M[:,1],V',M[:,3])
    Mz = hcat(M[:,1:2],V')
    result= @MVector [det(Mx), det(My), det(Mz)]
    result = result ./ nominator
end

function find_intersection(sl1::StraightLine,sl2::StraightLine)
    m = @MArray [sl1.step_vector[1] -1*sl2.step_vector[1];
    sl1.step_vector[2] -1*sl2.step_vector[2]]
    v = @MVector [sl2.reference_point[1] - sl1.reference_point[1],sl2.reference_point[2] - sl1.reference_point[2]]
    if iszero(det(m))
        m = @MArray [sl1.step_vector[1] -1*sl2.step_vector[1];
        sl1.step_vector[3] -1*sl2.step_vector[3]]
        v = @MVector [sl2.reference_point[1] - sl1.reference_point[1],sl2.reference_point[3] - sl1.reference_point[3]]
    else
        nothing
    end
    result = cramers_rule_2d(m,v)
    if isapprox(sl1.reference_point[3]+result[1]*sl1.step_vector[3],sl2.reference_point[3]+result[2]*sl2.step_vector[3])
        return [sl1.reference_point + result[1]*sl1.step_vector]
    else
        error("no intersection")
    end
end

function generate_random_startpositions(d::SolidStateDetector{T}, n::Int, Volume::NamedTuple=bounding_box(d), rng::AbstractRNG = MersenneTwister(), min_dist_from_boundary = 0.0001) where T
    delta = T(min_dist_from_boundary)
    n_filled::Int = 0
    positions = Vector{CartesianPoint{T}}(undef,n)
    while n_filled < n
        sample=CylindricalPoint{T}(rand(rng,Volume[:r_range].left:0.00001:Volume[:r_range].right),rand(rng,Volume[:φ_range].left:0.00001:Volume[:φ_range].right),rand(rng,Volume[:z_range].left:0.00001:Volume[:z_range].right))
        if !(sample in d.contacts) && contains(d,sample) && contains(d,CylindricalPoint{T}(sample.r+delta,sample.φ,sample.z))&& contains(d,CylindricalPoint{T}(sample.r-delta,sample.φ,sample.z))&& contains(d,CylindricalPoint{T}(sample.r,sample.φ,sample.z+delta))&& contains(d,CylindricalPoint{T}(sample.r,sample.φ,sample.z-delta))
            n_filled += 1
            positions[n_filled]=CartesianPoint(sample)
        end
    end
    positions
end

include("plot_recipes.jl")
