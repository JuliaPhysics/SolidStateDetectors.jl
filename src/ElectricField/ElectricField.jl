struct EField{T<:AbstractFloat}
    field::Array{SVector{3, T},3}

    function EField{T}(vectorfield) where T<:AbstractFloat
        return new{T}(vectorfield)
    end
end
function EField(vectorfield::Array{SVector{3,T},3}) where T<:AbstractFloat
    return EField{T}(vectorfield)
end



function get_magnitude_of_rθz_vector(vector::AbstractArray,cutoff=NaN)
    magn=0
    magn+=(vector[1])^2
    magn+=(vector[3])^2
    result = sqrt(magn)
    if isnan(cutoff)
        return result
    else
        if result >cutoff
            result=0
        end
        return result
    end
end


function get_electric_field_from_potential(ep::ElectricPotential{T}, fieldvector_coordinates=:xyz)::Array{SArray{Tuple{3},T,1,3}, 3} where {T <: AbstractFloat} 
    p = ep.data
    axr::Vector{T} = collect(ep.grid[:r])
    axθ::Vector{T} = collect(ep.grid[:θ])
    axz::Vector{T} = collect(ep.grid[:z])
    cyclic::T = ep.grid[:θ].interval.right
    ef = Array{SVector{3, T}}(undef, size(p)...)
    for iz in 1:size(ef, 3)
        for iθ in 1:size(ef, 2)
            for ir in 1:size(ef, 1)
                ### r ###
                if ir-1<1
                    Δp_r_1 = p[ir+1 ,iθ, iz] - p[ir ,iθ, iz]
                    d_r_1 = axr[ir+1]-axr[ir]
                    er = ( Δp_r_1 / d_r_1 )
                elseif  ir+1 > size(ef,1)
                    Δp_r_1 = p[ir ,iθ, iz]-p[ir-1 ,iθ, iz]
                    d_r_1 = axr[ir]-axr[ir-1]
                    er = ( Δp_r_1/d_r_1 )
                else
                    Δp_r_1 = p[ir+1 ,iθ, iz]-p[ir ,iθ, iz]
                    Δp_r_2 = p[ir ,iθ, iz]-p[ir-1 ,iθ, iz]
                    d_r_1 = axr[ir+1]-axr[ir]
                    d_r_2 = axr[ir]-axr[ir-1]
                    er = ( Δp_r_1/d_r_1 + Δp_r_2/d_r_2) / 2
                end
                ### θ ###
                if iθ < 2
                    Δp_θ_1 = p[ir ,iθ+1, iz]-p[ir ,iθ, iz]
                    Δp_θ_2 = p[ir ,iθ, iz]-p[ir ,end, iz]
                    d_θ_1 = (axθ[iθ+1]-axθ[iθ])/axr[ir]# to get the proper value in length units
                    d_θ_2 = (cyclic - axθ[end]) / axr[ir]
                    eθ = ( Δp_θ_1/d_θ_1 + Δp_θ_2/d_θ_2) / 2
                elseif iθ == size(ef,2)

                    Δp_θ_1 = p[ir ,1, iz]-p[ir ,iθ, iz]
                    Δp_θ_2 = p[ir ,iθ, iz]-p[ir ,iθ-1, iz]
                    d_θ_1 = (axθ[1]-axθ[iθ])/axr[ir]# to get the proper value in length units
                    d_θ_2 = (axθ[iθ]-axθ[iθ-1])/axr[ir]
                    eθ = ( Δp_θ_1/d_θ_1 + Δp_θ_2/d_θ_2) / 2
                else
                    Δp_θ_1 = p[ir ,iθ+1, iz]-p[ir ,iθ, iz]
                    Δp_θ_2 = p[ir ,iθ, iz]-p[ir ,iθ-1, iz]
                    d_θ_1 = (axθ[iθ+1]-axθ[iθ])/axr[ir]# to get the proper value in length units
                    d_θ_2 = (axθ[iθ]-axθ[iθ-1])/axr[ir]
                    eθ = ( Δp_θ_1/d_θ_1 + Δp_θ_2/d_θ_2) / 2
                end
                isinf(eθ) || isnan(eθ) ? eθ = 0.0 : nothing # for small radii and small distances(center of the grid) it would yield Infs or Nans
                if iz-1<1
                    Δp_z_1 = p[ir ,iθ, iz+1]-p[ir ,iθ, iz]
                    d_z_1 = axz[iz+1]-axz[iz]
                    ez = ( Δp_z_1/d_z_1 )
                elseif  iz+1 > size(ef,3)
                    Δp_z_1 = p[ir ,iθ, iz]-p[ir ,iθ, iz-1]
                    d_z_1 = axz[iz]-axz[iz-1]
                    ez = ( Δp_z_1/d_z_1 )
                else
                    Δp_z_1 = p[ir ,iθ, iz+1]-p[ir ,iθ, iz]
                    Δp_z_2 = p[ir ,iθ, iz]-p[ir ,iθ, iz-1]
                    d_z_1 = axz[iz+1]-axz[iz]
                    d_z_2 = axz[iz]-axz[iz-1]
                    ez = ( Δp_z_1/d_z_1 + Δp_z_2/d_z_2) / 2
                end
                e_vector = [-er,-eθ,-ez]
                ef[ir,iθ,iz] = e_vector
            end
        end
    end
    if fieldvector_coordinates == :xyz
        ef = convert_field_vectors_to_xyz(ef,axθ)
    end
    return ef
end
#
# function get_electric_field_from_potential(coax_grid::CylindricalGrid, point_type_array::PointTypes)
#     p=coax_grid.potential
#     ef = Array{Vector{Float32}}(undef,size(p,1),size(p,2),size(p,3))
#     for iz in 1:size(ef, 3)
#         for iθ in 1:size(ef, 2)
#             # isum = iz + iθ
#             for ir in 1:size(ef, 1)
#               #
#               # ninds = ir - 1, iθ - 1, iz - 1
#                       # rbinds = get_rb_inds(ninds...)
#                       # pwinds = rbinds .- 1
#               #
#               # pt = if iseven(isum + ir)
#               #   pw.point_types_even[pwinds...]
#               # else
#               #   pw.point_types_odd[pwinds...]
#               # end
#                 pt = point_type_array[ir, iθ ,iz]
#                 if pt & inside_bit == 0 ## means outside
#                     er=0.0
#                     dr=1.0
#                 elseif pt & r_lb_bit > 0 || ir-1 < 1
#                     dr = coax_grid.r[ir+1]-coax_grid.r[ir]
#                     er = p[ir+1 ,iθ, iz]-p[ir ,iθ, iz]
#                 elseif pt & r_rb_bit > 0 || ir+1 > size(ef,1)
#                     dr = coax_grid.r[ir]-coax_grid.r[ir-1]
#                     er = p[ir ,iθ, iz]-p[ir-1 ,iθ, iz]
#                 else
#                     dr = coax_grid.r[ir+1]-coax_grid.r[ir-1]
#                     er = p[ir+1 ,iθ, iz]-p[ir-1 ,iθ, iz]
#                 end
#
#                 #### θ ####
#                 if pt & inside_bit == 0 ## outside
#                     dθ=1.0
#                     eθ=0.0
#                 elseif iθ-1 < 1 ||pt & θ_lb_bit > 0
#                     dθ = coax_grid.θ[iθ+1]-coax_grid.θ[iθ]
#                     eθ = p[ir ,iθ+1, iz]-p[ir ,iθ, iz]
#                 elseif iθ+1 > size(ef,2) || pt & θ_rb_bit > 0
#                     dθ = coax_grid.θ[iθ]-coax_grid.θ[iθ-1]
#                     eθ = p[ir ,iθ, iz]-p[ir ,iθ-1, iz]
#                 else
#                     dθ = coax_grid.θ[iθ+1]-coax_grid.θ[iθ-1]
#                     eθ = p[ir ,iθ+1, iz]-p[ir ,iθ-1, iz]
#                 end
#
#                 #### z ####
#                 if pt & inside_bit == 0 ## outside
#                     dz=1.0
#                     ez=0.0
#                 elseif pt & z_lb_bit > 0 || iz-1 < 1## || (pt_z==fixed_point && pw.point_types[ 3, ir+1, iθ+1, iz]==outside)
#                     dz = coax_grid.z[iz+1]-coax_grid.z[iz]
#                     ez = p[ir ,iθ, iz+1]-p[ir ,iθ, iz]
#                     # @show ir, iθ, iz
#                     # @show ez, dz
#                 elseif pt & z_rb_bit > 0 || iz+1 > size(ef,3) ##|| (pt_z==fixed_point && pw.point_types[ 3, ir+1, iθ+1, iz+2]==outside)
#                     dz = coax_grid.z[iz]-coax_grid.z[iz-1]
#                     ez = p[ir ,iθ, iz]-p[ir ,iθ, iz-1]
#                 else
#                     dz = coax_grid.z[iz+1]-coax_grid.z[iz-1]
#                     ez = p[ir ,iθ, iz+1]-p[ir ,iθ, iz-1]
#                 end
#
#                 e_vector = [-er/dr,-eθ/dθ,-ez/dz]
#                 ef[ir,iθ,iz] = e_vector
#             end
#         end
#     end
#     # return ef
#     return EField(ef)
# end
function get_component_field(ef,component=:r,cutoff=NaN)
    components = [:r,:phi,:z]
    # component_index = findfirst(components,component)
    component_index = findfirst(x->x==component,components)
    ef_component = Array{Float32}(undef,size(ef,1),size(ef,2),size(ef,3))
    for iz in 1:size(ef, 3)
        for iθ in 1:size(ef, 2)
            for ir in 1:size(ef, 1)
                if !isnan(cutoff)
                    if abs(ef[ir,iθ ,iz][component_index]) >=cutoff
                        ef_component[ir,iθ,iz] = 0.0
                    else
                        ef_component[ir,iθ,iz] = ef[ir,iθ ,iz][component_index]
                    end
                else
                    ef_component[ir,iθ,iz] = ef[ir,iθ ,iz][component_index]
                end
            end
        end
    end
    return ef_component
end
function get_xyz_vector_from_rθz_vector(v::AbstractArray)::AbstractArray
    return [v[1]*cos(v[2]),v[1]*sin(v[2]),v[3]]
end

function get_xyz_vector_from_field_vector(field,r,θ,z,ir,iθ,iz)
    startpoint_vector = get_xyz_vector_from_rθz_vector([r,θ,z])
    endpoint_vector = get_xyz_vector_from_rθz_vector([r,θ,z]+field[ir,iθ,iz])
    xyz_vector = endpoint_vector-startpoint_vector
    for ic in 1:size(xyz_vector,1)
        isapprox(xyz_vector[ic],0.0) ? xyz_vector[ic] = 0.0 : nothing
    end
    return xyz_vector
end

function convert_field_vectors_to_xyz(field::Array{SArray{Tuple{3},T,1,3},3}, θa::Array{T, 1})::Array{SVector{3, T},3} where {T}
    field_xyz = Array{SVector{3,T},3}(undef, size(field)...);
    for (iθ, θ) in enumerate(θa)
        Rα::SMatrix{3,3,T} = @SArray([cos(θ) -1*sin(θ) 0;sin(θ) cos(θ) 0;0 0 1])
        for iz in axes(field, 3)
            for ir in axes(field, 1)
                field_xyz[ir,iθ,iz] = Rα * field[ir,iθ,iz]
                # field_xyz[ir,iθ,iz] = get_xyz_vector_from_field_vector(field[ir,iθ,iz], θ)
            end
        end
    end
    return field_xyz
end

function get_xyz_vector_from_field_vector(vector, α)
        # Rα = Array{AbstractFloat}(undef,2,2)
        # Rα[1,1]=cos(α)
        # Rα[1,2]=-1*sin(α)
        # Rα[2,1]=sin(α)
        # Rα[2,2]=cos(α)
        # result = Rα*vector[1:2]
        # push!(result,vector[3])
        # result

        Rα = @SArray([cos(α) -1*sin(α) 0;sin(α) cos(α) 0;0 0 1])
        result = Rα*vector
        result
end
function interpolated_scalarfield(ep::ScalarPotential{T}) where {T}
    knots = ep.grid.axes#(grid.r, grid.θ, grid.z)
    i = interpolate(knots, ep.data, Gridded(Linear()))
    vector_field_itp = extrapolate(i, Periodic())
    return vector_field_itp
end

function interpolated_vectorfield(vectorfield::AbstractArray{<:SVector{3, T},3}, ep::ElectricPotential{T}) where {T}
    knots = ep.grid.axes#(grid.r, grid.θ, grid.z)
    i = interpolate(knots, vectorfield, Gridded(Linear()))
    vectorfield_itp = extrapolate(i, Periodic())
    return vectorfield_itp
end

function setup_interpolated_efield(ef, ep::ScalarPotential{T}) where {T}# returns interpolation object
    knots = ep.grid.axes #(grid.r, grid.θ, grid.z)
    # itp = interpolate(knots,ef,Gridded(Linear()))
    itp = interpolate(knots, ef, Gridded(Linear()))
    return itp
end
function get_interpolated_efield_vector(itp, r, θ, z)
    return itp[r, θ, z]
end
function setup_interpolated_vectorfield(vectorfield, grid::CylindricalGrid{T}) where {T}
    knots = grid.axes #(grid.r, grid.θ, grid.z)
    i = interpolate(knots, vectorfield, Gridded(Linear()))
    vector_field_itp = extrapolate(i, Periodic())
    return vector_field_itp
end

function get_interpolated_drift_field(velocity_field, grid::CylindricalGrid{T}) where {T}
    knots = grid.axes #(grid.r, grid.θ, grid.z)
    # println(typeof(knots))
    # println(typeof(velocity_field))
    i=interpolate(knots, velocity_field, Gridded(Linear()))
    # println(typeof(i))
    velocity_field_itp = extrapolate(i, Periodic())
    return velocity_field_itp
end

include("plot_recipes.jl")
