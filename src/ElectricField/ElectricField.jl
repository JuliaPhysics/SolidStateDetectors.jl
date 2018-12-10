struct EField{T<:AbstractFloat}
    field::Array{SVector{3, T},3}

    function EField{T}(vectorfield) where T<:AbstractFloat
        return new{T}(vectorfield)
    end
end
function EField(vectorfield::Array{SVector{3,T},3}) where T<:AbstractFloat
    return EField{T}(vectorfield)
end


function get_magnitude_of_rφz_vector(vector::AbstractArray,cutoff=NaN)
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


function get_electric_field_from_potential(pg::CylindricalGrid{T},fieldvector_coordinates=:xyz)::Array{SArray{Tuple{3},T,1,3}, 3} where {T <: AbstractFloat} ## pg = potential_grid
    p=pg.potential
    ef = Array{SVector{3, T}}(undef, size(p)...)
    for iz in 1:size(ef, 3)
        for iφ in 1:size(ef, 2)
            for ir in 1:size(ef, 1)
                ### r ###
                if ir-1<1
                    Δp_r_1 = p[ir+1 ,iφ, iz] - p[ir ,iφ, iz]
                    d_r_1 = pg.r[ir+1]-pg.r[ir]
                    er = ( Δp_r_1 / d_r_1 )
                elseif  ir+1 > size(ef,1)
                    Δp_r_1 = p[ir ,iφ, iz]-p[ir-1 ,iφ, iz]
                    d_r_1 = pg.r[ir]-pg.r[ir-1]
                    er = ( Δp_r_1/d_r_1 )
                else
                    Δp_r_1 = p[ir+1 ,iφ, iz]-p[ir ,iφ, iz]
                    Δp_r_2 = p[ir ,iφ, iz]-p[ir-1 ,iφ, iz]
                    d_r_1 = pg.r[ir+1]-pg.r[ir]
                    d_r_2 = pg.r[ir]-pg.r[ir-1]
                    er = ( Δp_r_1/d_r_1 + Δp_r_2/d_r_2) / 2
                end
                ### φ ###
                if iφ < 2
                    Δp_φ_1 = p[ir ,iφ+1, iz]-p[ir ,iφ, iz]
                    Δp_φ_2 = p[ir ,iφ, iz]-p[ir ,end, iz]
                    d_φ_1 = (pg.φ[iφ+1]-pg.φ[iφ])/pg.r[ir]# to get the proper value in length units
                    d_φ_2 = (pg.cyclic - pg.φ[end]) / pg.r[ir]
                    eφ = ( Δp_φ_1/d_φ_1 + Δp_φ_2/d_φ_2) / 2
                elseif iφ == size(ef,2)

                    Δp_φ_1 = p[ir ,1, iz]-p[ir ,iφ, iz]
                    Δp_φ_2 = p[ir ,iφ, iz]-p[ir ,iφ-1, iz]
                    d_φ_1 = (pg.φ[1]-pg.φ[iφ])/pg.r[ir]# to get the proper value in length units
                    d_φ_2 = (pg.φ[iφ]-pg.φ[iφ-1])/pg.r[ir]
                    eφ = ( Δp_φ_1/d_φ_1 + Δp_φ_2/d_φ_2) / 2
                else
                    Δp_φ_1 = p[ir ,iφ+1, iz]-p[ir ,iφ, iz]
                    Δp_φ_2 = p[ir ,iφ, iz]-p[ir ,iφ-1, iz]
                    d_φ_1 = (pg.φ[iφ+1]-pg.φ[iφ])/pg.r[ir]# to get the proper value in length units
                    d_φ_2 = (pg.φ[iφ]-pg.φ[iφ-1])/pg.r[ir]
                    eφ = ( Δp_φ_1/d_φ_1 + Δp_φ_2/d_φ_2) / 2
                end
                isinf(eφ) || isnan(eφ) ? eφ = 0.0 : nothing # for small radii and small distances(center of the grid) it would yield Infs or Nans
                if iz-1<1
                    Δp_z_1 = p[ir ,iφ, iz+1]-p[ir ,iφ, iz]
                    d_z_1 = pg.z[iz+1]-pg.z[iz]
                    ez = ( Δp_z_1/d_z_1 )
                elseif  iz+1 > size(ef,3)
                    Δp_z_1 = p[ir ,iφ, iz]-p[ir ,iφ, iz-1]
                    d_z_1 = pg.z[iz]-pg.z[iz-1]
                    ez = ( Δp_z_1/d_z_1 )
                else
                    Δp_z_1 = p[ir ,iφ, iz+1]-p[ir ,iφ, iz]
                    Δp_z_2 = p[ir ,iφ, iz]-p[ir ,iφ, iz-1]
                    d_z_1 = pg.z[iz+1]-pg.z[iz]
                    d_z_2 = pg.z[iz]-pg.z[iz-1]
                    ez = ( Δp_z_1/d_z_1 + Δp_z_2/d_z_2) / 2
                end
                e_vector = [-er,-eφ,-ez]
                ef[ir,iφ,iz] = e_vector
            end
        end
    end
    if fieldvector_coordinates == :xyz
        ef = convert_field_vectors_to_xyz(ef,pg.φ)
    end
    return ef
end
#
# function get_electric_field_from_potential(coax_grid::CylindricalGrid, point_type_array::PointTypes)
#     p=coax_grid.potential
#     ef = Array{Vector{Float32}}(undef,size(p,1),size(p,2),size(p,3))
#     for iz in 1:size(ef, 3)
#         for iφ in 1:size(ef, 2)
#             # isum = iz + iφ
#             for ir in 1:size(ef, 1)
#               #
#               # ninds = ir - 1, iφ - 1, iz - 1
#                       # rbinds = get_rb_inds(ninds...)
#                       # pwinds = rbinds .- 1
#               #
#               # pt = if iseven(isum + ir)
#               #   pw.point_types_even[pwinds...]
#               # else
#               #   pw.point_types_odd[pwinds...]
#               # end
#                 pt = point_type_array[ir, iφ ,iz]
#                 if pt & inside_bit == 0 ## means outside
#                     er=0.0
#                     dr=1.0
#                 elseif pt & r_lb_bit > 0 || ir-1 < 1
#                     dr = coax_grid.r[ir+1]-coax_grid.r[ir]
#                     er = p[ir+1 ,iφ, iz]-p[ir ,iφ, iz]
#                 elseif pt & r_rb_bit > 0 || ir+1 > size(ef,1)
#                     dr = coax_grid.r[ir]-coax_grid.r[ir-1]
#                     er = p[ir ,iφ, iz]-p[ir-1 ,iφ, iz]
#                 else
#                     dr = coax_grid.r[ir+1]-coax_grid.r[ir-1]
#                     er = p[ir+1 ,iφ, iz]-p[ir-1 ,iφ, iz]
#                 end
#
#                 #### φ ####
#                 if pt & inside_bit == 0 ## outside
#                     dφ=1.0
#                     eφ=0.0
#                 elseif iφ-1 < 1 ||pt & φ_lb_bit > 0
#                     dφ = coax_grid.φ[iφ+1]-coax_grid.φ[iφ]
#                     eφ = p[ir ,iφ+1, iz]-p[ir ,iφ, iz]
#                 elseif iφ+1 > size(ef,2) || pt & φ_rb_bit > 0
#                     dφ = coax_grid.φ[iφ]-coax_grid.φ[iφ-1]
#                     eφ = p[ir ,iφ, iz]-p[ir ,iφ-1, iz]
#                 else
#                     dφ = coax_grid.φ[iφ+1]-coax_grid.φ[iφ-1]
#                     eφ = p[ir ,iφ+1, iz]-p[ir ,iφ-1, iz]
#                 end
#
#                 #### z ####
#                 if pt & inside_bit == 0 ## outside
#                     dz=1.0
#                     ez=0.0
#                 elseif pt & z_lb_bit > 0 || iz-1 < 1## || (pt_z==fixed_point && pw.point_types[ 3, ir+1, iφ+1, iz]==outside)
#                     dz = coax_grid.z[iz+1]-coax_grid.z[iz]
#                     ez = p[ir ,iφ, iz+1]-p[ir ,iφ, iz]
#                     # @show ir, iφ, iz
#                     # @show ez, dz
#                 elseif pt & z_rb_bit > 0 || iz+1 > size(ef,3) ##|| (pt_z==fixed_point && pw.point_types[ 3, ir+1, iφ+1, iz+2]==outside)
#                     dz = coax_grid.z[iz]-coax_grid.z[iz-1]
#                     ez = p[ir ,iφ, iz]-p[ir ,iφ, iz-1]
#                 else
#                     dz = coax_grid.z[iz+1]-coax_grid.z[iz-1]
#                     ez = p[ir ,iφ, iz+1]-p[ir ,iφ, iz-1]
#                 end
#
#                 e_vector = [-er/dr,-eφ/dφ,-ez/dz]
#                 ef[ir,iφ,iz] = e_vector
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
                for iφ in 1:size(ef, 2)
                        for ir in 1:size(ef, 1)
                                if !isnan(cutoff)
                                        if abs(ef[ir,iφ ,iz][component_index]) >=cutoff
                                                ef_component[ir,iφ,iz] = 0.0
                                        else
                                                ef_component[ir,iφ,iz] = ef[ir,iφ ,iz][component_index]
                                        end
                                else
                                        ef_component[ir,iφ,iz] = ef[ir,iφ ,iz][component_index]
                                end
                        end
                end
        end
        return ef_component
end
function get_xyz_vector_from_rϕz_vector(v::AbstractArray)::AbstractArray
        return [v[1]*cos(v[2]),v[1]*sin(v[2]),v[3]]
end
function get_xyz_vector_from_rφz_vector(v::AbstractArray)::AbstractArray
        [v[1]*cos(v[2]),v[1]*sin(v[2]),v[3]]
end

function get_xyz_vector_from_field_vector(field,r,ϕ,z,ir,iϕ,iz)
        startpoint_vector = get_xyz_vector_from_rϕz_vector([r,ϕ,z])
        endpoint_vector = get_xyz_vector_from_rϕz_vector([r,ϕ,z]+field[ir,iϕ,iz])
        xyz_vector = endpoint_vector-startpoint_vector
        for ic in 1:size(xyz_vector,1)
                isapprox(xyz_vector[ic],0.0) ? xyz_vector[ic] = 0.0 : nothing
        end
        return xyz_vector
end

function convert_field_vectors_to_xyz(field::Array{SArray{Tuple{3},T,1,3},3}, φa::Array{T, 1})::Array{SVector{3, T},3} where {T}
        field_xyz = Array{SVector{3,T},3}(undef, size(field)...);
        for (iφ, φ) in enumerate(φa)
                Rα::SMatrix{3,3,T} = @SArray([cos(φ) -1*sin(φ) 0;sin(φ) cos(φ) 0;0 0 1])
                for iz in axes(field, 3)
                        for ir in axes(field, 1)
                                field_xyz[ir,iφ,iz] = Rα * field[ir,iφ,iz]
                                # field_xyz[ir,iφ,iz] = get_xyz_vector_from_field_vector(field[ir,iφ,iz], φ)
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
function interpolated_scalarfield(grid::CylindricalGrid)
    knots = (grid.r,grid.φ,grid.z)
    i=interpolate(knots,grid.potential,Gridded(Linear()))
    vector_field_itp =extrapolate(i,Periodic())
    return vector_field_itp
end

function interpolated_vectorfield(vectorfield::AbstractArray{<:SVector{3,T},3},grid::CylindricalGrid{T}) where T<:Real
    knots = (grid.r,grid.φ,grid.z)
    i=interpolate(knots,vectorfield,Gridded(Linear()))
    vectorfield_itp =extrapolate(i,Periodic())
    return vectorfield_itp
end

function setup_interpolated_efield(ef, grid::CylindricalGrid) # returns interpolation object
    knots = (grid.r,grid.φ,grid.z)
    # itp = interpolate(knots,ef,Gridded(Linear()))
    itp = interpolate(knots,ef,Gridded(Linear()))
    return itp
end
function get_interpolated_efield_vector(itp,r,φ,z)
    return itp[r,φ,z]
end
function setup_interpolated_vectorfield(vectorfield, grid::CylindricalGrid)
    knots = (grid.r,grid.φ,grid.z)
    i=interpolate(knots,vectorfield,Gridded(Linear()))
    vector_field_itp =extrapolate(i,Periodic())
    return vector_field_itp
end

function get_interpolated_drift_field(velocity_field, grid::CylindricalGrid)
    knots = (grid.r,grid.φ,grid.z)
    # println(typeof(knots))
    # println(typeof(velocity_field))
    i=interpolate(knots,velocity_field,Gridded(Linear()))
    # println(typeof(i))
    velocity_field_itp =extrapolate(i,Periodic())
    return velocity_field_itp
end
include("plot_recipes.jl")
