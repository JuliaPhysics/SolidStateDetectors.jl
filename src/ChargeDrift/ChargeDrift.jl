
# function drift_charge!(drift_path::AbstractVector{<:SVector{3,<:Real}},det::SolidStateDetector, startpos::SVector{3,<:Real}, delta_t::Real, velocity_field::AbstractArray{<:SVector{3,<:Real},3})::Real
function drift_charge!(
    drift_path::AbstractVector{<:SVector{3,T}},
    det::SolidStateDetector,
    startpos::SVector{3,<:Real},
    delta_t::Real,
    velocity_field::Interpolations.Extrapolation{<:SVector{3,<:Real},3}
) where T <:Real
    drifttime::T = 0.0
    done::Bool = false
    stepvector = @SVector zeros(T,3)
    pos_cyl = Cylindrical{T}(0,0,0)
    crossing_pos_cyl = Cylindrical{T}(0,0,0)
    crossing_pos_type::Symbol = :Idle
    boundary_index::Int=0

    drift_path[1] = startpos
    for istep in eachindex(drift_path)[2:end]
        if done == false
            pos_cyl = CylindricalPoint(drift_path[istep-1]) # update pos_cyl
            if contains(det, pos_cyl)
                stepvector = getvelocityvector(velocity_field, pos_cyl) * delta_t
                drift_path[istep] = drift_path[istep-1] + stepvector
            elseif istep < 3
                done = true
                drifttime = delta_t * (istep-1)
                throw(ErrorException("start position $pos_cyl is not contained in detector"))
            else
                crossing_pos_cyl, crossing_pos_type, boundary_index = get_crossing_pos(det, drift_path[istep-2], drift_path[istep-1])
                if crossing_pos_type == :electrode
                    done = true
                    drifttime = delta_t * (istep-1)
                    drift_path[istep-1] = CartFromCyl(crossing_pos_cyl)
                    drift_path[istep] = drift_path[istep-1]
                elseif crossing_pos_type == :floating_boundary
                    stepvector = getvelocityvector(velocity_field, crossing_pos_cyl) * delta_t
                    projected_vector = project_to_plane(stepvector, get_fb_plane(det, crossing_pos_cyl, boundary_index).n⃗)
                    drift_path[istep-1] = CartFromCyl(crossing_pos_cyl)
                    drift_path[istep] = drift_path[istep-1] + projected_vector
                else
                    @show crossing_pos_cyl,crossing_pos_type, boundary_index
                    @show typeof(crossing_pos_cyl)
                    @show typeof(pos_cyl)
                    throw(ErrorException("Internal error"))
                end
            end
        elseif done == true
            drift_path[istep] = drift_path[istep-1]
        end
    end
    drifttime
end


function drift_and_acc_charge!(
    drift_path::AbstractVector{<:SomeCartesianPoint{<:Real}},
    charge_signal::AbstractVector{<:RealQuantity},
    detector::SolidStateDetector,
    velocity_field::Interpolations.Extrapolation{<:SVector{3},3},
    weighting_potential::Interpolations.Extrapolation{<:Real,3},
    startpos::SomeCartesianPoint{<:RealQuantity},
    charge::RealQuantity,
    delta_t::RealQuantity,
)
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


# function CylFromCart!(p_cyl::Cylindrical{T}, p::SVector{3,T})::Cylindrical where T <:Real
#     t=CylindricalFromCartesian()(p)
#     θ=t.θ
#     while θ<0
#         θ+=2π
#     end
#     while θ>2π
#         θ-=2π
#     end
#     p_cyl.r= t.r
#     p_cyl.θ = θ
#     p_cyl.z = t.z
#     nothing
# end

function driftonecharge(det::SolidStateDetector{T}, startposition::AbstractArray, nsteps::Int, steptime, charge::Symbol, interpolation_field_e, interpolation_field_h) where {T} #where {T<:AbstractFloat}#{T=get_precision_type(det)}
    # T::Type = eltype(startposition)
    pathx::AbstractArray{T, 1} = zeros(T, nsteps)
    pathy::AbstractArray{T, 1} = zeros(T, nsteps)
    pathz::AbstractArray{T, 1} = zeros(T, nsteps)

    pathx[1] = startposition[1]
    pathy[1] = startposition[2]
    pathz[1] = startposition[3]
    prevposition::SVector{3, T} = startposition
    drifttime::T = 0.0
    done::Bool = false
    for istep in 2:nsteps
        if done == false
            stepvector = zeros(T, 3)
            position_cyl::Cylindrical{T}=CylFromCart(prevposition)
            if contains(det, position_cyl)
                if charge == :e
                    stepvector = getvelocityvector(interpolation_field_e,prevposition) * steptime
                elseif charge == :h
                    stepvector = getvelocityvector(interpolation_field_h,prevposition) * steptime
                else
                    error("wrong symbol for charge")
                end
                nextposition = prevposition + stepvector
                pathx[istep] = nextposition[1]
                pathy[istep] = nextposition[2]
                pathz[istep] = nextposition[3]
                prevposition = nextposition
            elseif istep < 3
                done = true
                drifttime = istep-1*steptime
                @warn "spawn_position is not contained $position_cyl"
                pathx[istep] = pathx[istep-1]
                pathy[istep] = pathy[istep-1]
                pathz[istep] = pathz[istep-1]
            else
                crossingpoint_cyl, cp_type, boundary_index = get_cyl_intersection_qad2(det,SVector{3,T}(pathx[istep-2],pathy[istep-2],pathz[istep-2]),prevposition)
                if cp_type == :electrode
                    done = true
                    drifttime = istep-1*steptime
                    crossingpoint=CartFromCyl(crossingpoint_cyl)
                    pathx[istep-1] = crossingpoint[1]
                    pathy[istep-1] = crossingpoint[2]
                    pathz[istep-1] = crossingpoint[3]

                    pathx[istep] = crossingpoint[1]
                    pathy[istep] = crossingpoint[2]
                    pathz[istep] = crossingpoint[3]

                elseif cp_type == :floating_boundary
                    println("Boundary Handling triggered...")
                    if charge == :e
                        stepvector = getvelocityvector(interpolation_field_e,prevposition) * steptime
                    elseif charge == :h
                        stepvector = getvelocityvector(interpolation_field_h,prevposition) * steptime
                    else
                        error("wrong symbol for charge")
                    end
                    projected_vector::SVector{3,T} = project_to_plane(stepvector,get_fb_plane(det,crossingpoint_cyl,boundary_index).n⃗)
                    crossingpoint=CartFromCyl(crossingpoint_cyl)
                    pathx[istep-1] = crossingpoint[1]
                    pathy[istep-1] = crossingpoint[2]
                    pathz[istep-1] = crossingpoint[3]

                    nextposition = crossingpoint + projected_vector
                    pathx[istep] = nextposition[1]
                    pathy[istep] = nextposition[2]
                    pathz[istep] = nextposition[3]
                    prevposition = nextposition
                else
                    @show crossingpoint_cyl,cp_type, boundary_index
                    @show typeof(crossingpoint_cyl)
                    @show position_cyl
                    @show typeof(position_cyl)
                    error("sth went wrong: ...")
                end
            end
        elseif done == true
            pathx[istep] = pathx[istep - 1]
            pathy[istep] = pathy[istep - 1]
            pathz[istep] = pathz[istep - 1]
        end
    end

    return pathx, pathy, pathz, drifttime
end

function driftonecharge_wo(det, startposition::AbstractArray, nsteps::Integer, steptime::AbstractFloat, charge::Symbol, interpolation_field_e, interpolation_field_h)
    # T::Type = eltype(startposition)
    T::Type = eltype(det.crystal_radius)
    pathx = zeros(T, nsteps)
    pathy = zeros(T, nsteps)
    pathz = zeros(T, nsteps)

    pathx[1] = startposition[1]
    pathy[1] = startposition[2]
    pathz[1] = startposition[3]
    prevposition = startposition
    # println("test")
    done::Bool=false
    for istep in 2:nsteps
        # println("$istep,$done")
        if done==false
            stepvector = zeros(T, 3)
            position_cyl=CylFromCart(prevposition)
            if contains(det, position_cyl)
                # println("in$istep")
            # rθz_position=get_rθz_from_xyz(prevposition)
            # if contains(det, rθz_position...)# contains in rθz position returns true if inside detector

            # if contains(det, (prevposition[1], prevposition[2], prevposition[3])) == 1
                if charge == :e
                    stepvector = getvelocityvector(interpolation_field_e,prevposition) * steptime
                elseif charge == :h
                    stepvector = getvelocityvector(interpolation_field_h,prevposition) * steptime
                else
                    error("wrong symbol for charge")
                end
            end
            nextposition = prevposition + stepvector
            pathx[istep] = nextposition[1]
            pathy[istep] = nextposition[2]
            pathz[istep] = nextposition[3]
            prevposition = nextposition

        elseif done == true
            pathx[istep] = pathx[istep-1]
            pathy[istep] = pathy[istep-1]
            pathz[istep] = pathz[istep-1]
        end
    end

    return pathx, pathy, pathz
end
function driftonecharge(det::SolidStateDetector, startposition_cyl::Cylindrical, nsteps::Integer, steptime::AbstractFloat, charge::Symbol, interpolation_field_e, interpolation_field_h)
    T=typeof(det.crystal_radius)
    driftonecharge(det, CartFromCyl(Cylindrical{T}(startposition_cyl.r,startposition_cyl.θ,startposition_cyl.z)),
                nsteps, T(steptime), charge, interpolation_field_e, interpolation_field_h)
end

struct DriftPath{T<:Real}
    e_path::AbstractVector{SVector{3,T}}
    h_path::AbstractVector{SVector{3,T}}

    DriftPath{T}(e_path::AbstractVector{<:SVector{3,T}},h_path::AbstractVector{<:SVector{3,T}}) where {T <:Real} = new{T}(e_path,h_path)
end

function DriftPath(e_path::AbstractVector{<:SVector{3,T}},h_path::AbstractVector{<:SVector{3,T}}) where T <:Real
    return DriftPath{T}(e_path,h_path)
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

function drift_charges(detector::SolidStateDetector, starting_positions::AbstractArray{SVector{3,T}}, velocity_field_e::Interpolations.Extrapolation{SVector{3,Float64},3}, velocity_field_h::Interpolations.Extrapolation{SVector{3,Float64},3}; delta_t::T=T(1f-9), n_steps::Int = 2000) where T <: Real
    drift_paths  = Vector{DriftPath{T}}(undef,size(starting_positions,1))
    drift_path_e = [@SVector zeros(T,3) for i in 1:n_steps]
    drift_path_h = [@SVector zeros(T,3) for i in 1:n_steps]

    for j in eachindex(starting_positions)
        drift_charge!(drift_path_e, detector, starting_positions[j], delta_t, velocity_field_e)
        drift_charge!(drift_path_h, detector, starting_positions[j], delta_t, velocity_field_h)
        drift_paths[j] = DriftPath(deepcopy(drift_path_e), deepcopy(drift_path_h))
    end

    drift_paths
end
function drift_charges(detector::SolidStateDetector, starting_positions::AbstractArray{SVector{3,T}}, velocity_field_e::Interpolations.Extrapolation{SVector{3,Float64},3}, velocity_field_h::Interpolations.Extrapolation{SVector{3,Float64},3},energy_depositions::AbstractVector{T},weighting_potentials::AbstractVector{<:Interpolations.Extrapolation{T,3}}; delta_t::T=T(1f-9), n_steps::Int = 2000) where T <: Real
    drift_paths =  drift_charges(detector, starting_positions, velocity_field_e, velocity_field_h)
    signal=Vector{AbstractVector{T}}(undef,size(weighting_potentials,1))
    for i in eachindex(weighting_potentials)
        signal[i]=pulse_from_drift_paths(drift_paths, energy_depositions, weighting_potentials[i])
    end
    return ChargeDriftEvent(drift_paths,signal)
end

function get_crossing_pos( detector::SolidStateDetector,point_in::SVector{3,T},point_out::SVector{3,T}; dig::Int=6, max_n_iter::Int=2000) where T <:Real
    T==Float64 ? dig = 13 : nothing
    current_pos_in = MVector{3,T}(point_in...)
    current_pos_out = MVector{3,T}(point_out...)
    current_pos_xyz::SVector{3,T} = 0.5*(current_pos_in .+ current_pos_out)
    current_pos_cyl::Cylindrical{T} = CylFromCart(current_pos_xyz)
    i_limit::Int = 0
    while (point_type(detector, current_pos_cyl.r, current_pos_cyl.θ, current_pos_cyl.z)[2]==0) && i_limit < max_n_iter
        if contains(detector,current_pos_cyl)
            current_pos_in[1] = current_pos_xyz[1]
            current_pos_in[2] = current_pos_xyz[2]
            current_pos_in[3] = current_pos_xyz[3]
        else
            current_pos_out[1] = current_pos_xyz[1]
            current_pos_out[2] = current_pos_xyz[2]
            current_pos_out[3] = current_pos_xyz[3]
        end
        current_pos_xyz = 0.5*(current_pos_out .+ current_pos_in)
        current_pos_cyl=CylFromCart(current_pos_xyz)
        i_limit += 1
    end
    crossing_pos_type, boundary_index = point_type(detector,current_pos_cyl.r,current_pos_cyl.θ,current_pos_cyl.z)
    if i_limit==max_n_iter
        # println("Boundary Point was rounded.")
        # @show round(current_pos_cyl.r,digits=dig),round(current_pos_cyl.θ,digits=dig),round(current_pos_cyl.z,digits=dig)
        crossing_pos_type, boundary_index = point_type(detector,round(current_pos_cyl.r,digits=dig),round(current_pos_cyl.θ,digits=dig),round(current_pos_cyl.z,digits=dig))
        # @show crossing_pos_type, boundary_index
    end
    current_pos_cyl, crossing_pos_type, boundary_index
end

function get_cyl_intersection_qad2(det::SolidStateDetector, p_in::AbstractArray, p_out::AbstractArray, dig=6, max_n_iter::Int=2000)
    T = typeof(det.crystal_length)
    T == Float64 ? dig=11 : nothing
    pos::SVector{3,T} = 0.5*(p_out .+ p_in)
    pos_cyl::Cylindrical{T}=CylFromCart(pos)
    i_limit=0
    while (point_type(det,pos_cyl.r,pos_cyl.θ,pos_cyl.z)[2]==0) && i_limit < max_n_iter
        contains(det,pos_cyl) ? p_in=pos : p_out = pos
        pos = 0.5*(p_out .+ p_in)
        pos_cyl=CylFromCart(pos)
        i_limit+=1
    end

    mypointtype, index = point_type(det,pos_cyl.r,pos_cyl.θ,pos_cyl.z)
    if i_limit==max_n_iter
        println("Boundary Point was rounded.")
        @show round(pos_cyl.r,digits=dig),round(pos_cyl.θ,digits=dig),round(pos_cyl.z,digits=dig)
        mypointtype, index = point_type(det,round(pos_cyl.r,digits=dig),round(pos_cyl.θ,digits=dig),round(pos_cyl.z,digits=dig))
        @show mypointtype, index
    end
    return Cylindrical{T}(round(pos_cyl.r,digits=dig),round(pos_cyl.θ,digits=dig),round(pos_cyl.z,digits=dig)), mypointtype, index
    # Alternative Outputs, predefined
    # return Cylindrical(pos_cyl.r,pos_cyl.θ,pos_cyl.z)
    # return SVector(round(pos[1],sigdigits=sigdig),round(pos[2],sigdigits=sigdig),round(pos[3],sigdigits=sigdig))
end

function project_to_plane(v⃗::AbstractArray,n⃗::AbstractArray) #Vector to be projected, #normal vector of plane
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

function CartFromCyl(p::Cylindrical)::AbstractArray
    CartesianFromCylindrical()(p)
end

# TODO: Improve this (slow and hacky)
function CylFromCart(p::SVector{3,T})::Cylindrical where T <:AbstractFloat
    t=CylindricalFromCartesian()(p)
    θ=t.θ
    while θ<0
        θ+=2π
    end
    while θ>2π
        θ-=2π
    end
    if T==Float32
        return Cylindrical{T}(round(t.r,digits=7),θ,round(t.z,digits=7))
    else
        return Cylindrical{T}(t.r,θ,t.z)
    end
end


function get_cartesian_paths(a::Array{CylindricalPoint,1})
    xpath = Array{eltype(a[1].r),1}(undef,size(a))
    ypath = Array{eltype(a[1].r),1}(undef,size(a))
    zpath = Array{eltype(a[1].r),1}(undef,size(a))
    xyz = map(x->CartesianPoint(x),a)
    for i in eachindex(xyz)
        xpath[i]=xyz[i].x
        ypath[i]=xyz[i].y
        zpath[i]=xyz[i].z
    end
    xpath,ypath,zpath
end
function get_cartesian_paths(a::Array{SArray{Tuple{3},T,1,3},1};scaling::Real = 1) where T <:Real
    xpath = Array{T,1}(undef,size(a))
    ypath = Array{T,1}(undef,size(a))
    zpath = Array{T,1}(undef,size(a))
    for i in eachindex(a)
        xpath[i]=a[i][1]*scaling
        ypath[i]=a[i][2]*scaling
        zpath[i]=a[i][3]*scaling
    end
    xpath,ypath,zpath
end
function get_rθz_vector_from_xyz(vector::AbstractArray)::AbstractArray
  @fastmath begin
    x = vector[1]
    y = vector[2]
    z = vector[3]
    r = sqrt(x^2 + y^2)
    # θ = atan(y , x)
    θ = atan(y / x)
    x <  0 && y >= 0 ? θ += π  : nothing
    x <= 0 && y <  0 ? θ += π  : nothing
    x >  0 && y <  0 ? θ += 2π : nothing
  end
  return [r, θ, z]
end

function get_rθz_from_xyz(v::AbstractArray)::AbstractArray
    r = sqrt(v[1]^2+v[2]^2)
    θ = atan(v[2],v[1])
    while θ<0
        θ+=2π
    end
    while θ>2π
        θ-=2π
    end
    z=v[3]
    [r,θ,z]
end

function convert_xyz_vectorfield_in_rθz_to_rθz_vectorfield_in_rθz(field,grid)
    outputfield = Array{Array{Float32,1}}(undef,size(field,1),size(field,2),size(field,3))
    for (iz,z) in enumerate(grid.z)
        for (iθ,θ) in enumerate(grid.θ)
            for (ir,r) in enumerate(grid.r)
                outputfield[ir,iθ,iz]=get_rθz_from_xyz(field[ir,iθ,iz])
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
    θ = v1[2]+atan(v2[1]*sin(v2[2]-v1[2]),v1[1]+v2[1]*cos(v2[2]-v1[2]))
    z = v1[3]+v2[3]
    [r,θ,z]
end


@inline getvelocityvector(interpolation_field, position::SomeCartesianPoint) =
    getvelocityvector(interpolation_field, convert(CylindricalPoint, position))

@inline getvelocityvector(interpolated_vectorfield, point::Cylindrical) =
    getvelocityvector(interpolated_vectorfield, CylindricalPoint(point))

@inline getvelocityvector(interpolated_vectorfield, point::CylindricalPoint{T}) where T =
    interpolated_vectorfield(point.r, point.θ, point.z)::SVector{3,T}


struct Trajectory
    startposition::Vector
    n_charges::Int
    indv_charge::Int
    n_steps::Int
    timestep::AbstractFloat
    rθz_path::Array{Cylindrical,1}
    x_path::Vector
    y_path::Vector
    z_path::Vector
    endposition::Vector
    drifttime::AbstractFloat

    function Trajectory(startposition, n_charges, indv_charge, n_steps, timestep, rθz_path, x_path, y_path, z_path,endposition, drifttime)
        return new(startposition, n_charges, indv_charge, n_steps, timestep, rθz_path, x_path,y_path, z_path, endposition, drifttime)
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

function Trajectory(det::SolidStateDetector,startposition_rθz::Cylindrical,n_steps::Int,timestep,carriertype::Symbol,velocity_field,n_charges = 1)
    x_path,  y_path, z_path, drifttime = driftonecharge(det,startposition_rθz,n_steps,timestep,carriertype,velocity_field,velocity_field)
    rθz_path =[]
    carriertype == :e ? charge = -1 : charge = 1
    Trajectory(CartFromCyl(startposition_rθz), n_charges, charge, n_steps, timestep, rθz_path, x_path, y_path, z_path, [x_path[end],y_path[end],z_path[end]], drifttime)
end

function Trajectory(det::SolidStateDetector,startposition::AbstractArray,n_steps::Int,timestep,carriertype::Symbol,velocity_field, n_charges = 1)
    x_path,  y_path, z_path, drifttime = driftonecharge(det,startposition,n_steps,timestep,carriertype,velocity_field,velocity_field)
    rθz_path =[]
    carriertype == :e ? charge = -1 : charge = 1
    Trajectory(startposition, n_charges, charge, n_steps, timestep, rθz_path, x_path, y_path, z_path, [x_path[end],y_path[end],z_path[end]], drifttime)
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
    pos_cyl::Cylindrical{T} = CylFromCart(pos)
    # isapprox(Wpot_interp(pos_cyl.r,pos_cyl.θ,pos_cyl.z),0.0,atol = 5e-5) ? T(0.0) : Wpot_interp(pos_cyl.r,pos_cyl.θ,pos_cyl.z)
    Wpot_interp(pos_cyl.r,pos_cyl.θ,pos_cyl.z)
end


function accumulate_charge!(charge_signal::AbstractVector{<:RealQuantity}, drift_path::AbstractVector{<:SVector{3,T}}, Wpot_interp::Interpolations.Extrapolation{T,3}, charge::RealQuantity) where T <:Real
    idxs = eachindex(charge_signal)
    eachindex(drift_path) == idxs || throw(ArgumentError("Vectors charge_signal and drift_path must have same indices"))
    pos_prev = drift_path[first(idxs)]
    V_0 = _get_wpot_at(Wpot_interp, pos_prev)
    V_prev::typeof(V_0) = V_0
    @inbounds for i in idxs
        pos_current = drift_path[i]
        V_current = pos_current ≈ pos_prev ? V_prev : _get_wpot_at(Wpot_interp, drift_path[i])
        charge_signal[i] = muladd(charge, V_current - V_0, charge_signal[i])
        V_prev = V_current
        pos_prev = pos_current
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
    hit_pos::AbstractVector{<:AbstractVector{<:SomeCartesianPoint}},
    hit_edep::AbstractVector{<:AbstractVector{<:RealQuantity}},
    E_ionisation::RealQuantity,
    delta_t::RealQuantity,
)
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
            unitless_startpos = to_internal_units(u"m", startpos_vec[i_hit])
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


function Pulse(Trajectory_e::Trajectory, Trajectory_h::Trajectory, wp_int, energy::T = 1.0) where T <: AbstractFloat#wp_int = interpolation field for the weighting potential
    n_steps = Trajectory_e.n_steps
    electron_contribution::Array{T,1} = zeros(T, n_steps)
    hole_contribution::Array{T,1} = zeros(T, n_steps)
    signal::Array{T,1} = zeros(T,n_steps)
    for i_step in eachindex(electron_contribution)
        curr_loc_xyz_e::SVector{3,T} = SVector(Trajectory_e.x_path[i_step],Trajectory_e.y_path[i_step],Trajectory_e.z_path[i_step])
        # curr_loc_rθz_e=get_rθz_from_xyz(curr_loc_xyz_e)
        curr_loc_rθz_e::Cylindrical = CylFromCart(curr_loc_xyz_e)
        curr_loc_xyz_h::SVector{3,T} = SVector(Trajectory_h.x_path[i_step],Trajectory_h.y_path[i_step],Trajectory_h.z_path[i_step])
        # curr_loc_rθz_h=get_rθz_from_xyz(curr_loc_xyz_h)
        curr_loc_rθz_h::Cylindrical=CylFromCart(curr_loc_xyz_h)
        # electron_contribution[i_step] = -1 * energy * wp_int(curr_loc_rθz_e...)
        electron_contribution[i_step] = -1 * energy * wp_int(curr_loc_rθz_e.r, curr_loc_rθz_e.θ, curr_loc_rθz_e.z)
        hole_contribution[i_step] = energy * wp_int(curr_loc_rθz_h.r, curr_loc_rθz_h.θ, curr_loc_rθz_h.z)
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
    function Energy_Deposition(location::Cylindrical,energy::AbstractFloat,n_eh_pairs::Int)
        return new(CartFromCyl(location),energy,n_eh_pairs)
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

# struct ChargeDriftEvent
#     drift_paths_e
#     drift_paths_h
#     pulses
# end
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
function get_fb_plane(det,pCyl::Cylindrical,boundary_index)::Plane
    T=typeof(det.crystal_length)
    pCart::SVector{3,T}=CartFromCyl(pCyl)
    if det.floating_boundary_r_ranges[boundary_index][1]==det.floating_boundary_r_ranges[boundary_index][2]
        # println("radial boundary")
        return Plane(pCart,SVector{3,T}(get_xyz_vector_from_field_vector([1,0,0],pCyl.θ)...))
    elseif det.floating_boundary_phi_ranges[boundary_index][1]==det.floating_boundary_phi_ranges[boundary_index][2]
        # println("phi boundary")
        return Plane(pCart,SVector{3,T}(get_xyz_vector_from_field_vector([0,1,0],pCyl.θ+π/2)...))
    elseif det.floating_boundary_z_ranges[boundary_index][1]==det.floating_boundary_z_ranges[boundary_index][2]
        # println("z boundary")
        return Plane(pCart, SVector{3,T}(0,0,1))
    end
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

include("plot_recipes.jl")
