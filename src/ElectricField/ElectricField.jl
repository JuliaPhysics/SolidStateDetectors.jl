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


function get_electric_field_from_potential(ep::ElectricPotential{T, 3, :Cylindrical}, pointtypes::PointTypes{T}, fieldvector_coordinates=:xyz)::Array{SArray{Tuple{3},T,1,3}, 3} where {T <: AbstractFloat}
    p = ep.data
    axr::Vector{T} = collect(ep.grid[:r])
    axφ::Vector{T} = collect(ep.grid[:φ])
    axz::Vector{T} = collect(ep.grid[:z])

    cyclic::T = ep.grid[:φ].interval.right
    ef = Array{SVector{3, T}}(undef, size(p)...)
    for iz in 1:size(ef, 3)
        for iφ in 1:size(ef, 2)
            for ir in 1:size(ef, 1)
                ### r ###
                if ir-1<1
                    Δp_r_1 = p[ir+1 ,iφ, iz] - p[ir ,iφ, iz]
                    d_r_1 = axr[ir+1]-axr[ir]
                    er = ( Δp_r_1 / d_r_1 )
                elseif  ir+1 > size(ef,1)
                    Δp_r_1 = p[ir ,iφ, iz]-p[ir-1 ,iφ, iz]
                    d_r_1 = axr[ir]-axr[ir-1]
                    er = ( Δp_r_1/d_r_1 )
                else
                    Δp_r_1 = p[ir+1 ,iφ, iz]-p[ir ,iφ, iz]
                    Δp_r_2 = p[ir ,iφ, iz]-p[ir-1 ,iφ, iz]
                    d_r_1 = axr[ir+1]-axr[ir]
                    d_r_2 = axr[ir]-axr[ir-1]
                    er = ( Δp_r_1/d_r_1 + Δp_r_2/d_r_2) / 2
                end
                ### φ ###
                if iφ < 2
                    Δp_φ_1 = p[ir ,iφ+1, iz]-p[ir ,iφ, iz]
                    Δp_φ_2 = p[ir ,iφ, iz]-p[ir ,end, iz]
                    d_φ_1 = (axφ[iφ+1]-axφ[iφ]) * axr[ir]# to get the proper value in length units
                    d_φ_2 = (cyclic - axφ[end]) * axr[ir]
                    eφ = ( Δp_φ_1/d_φ_1 + Δp_φ_2/d_φ_2) / 2
                elseif iφ == size(ef,2)
                    Δp_φ_1 = p[ir ,1, iz]-p[ir ,iφ, iz]
                    Δp_φ_2 = p[ir ,iφ, iz]-p[ir ,iφ-1, iz]
                    d_φ_1 = (axφ[1]-axφ[iφ]) * axr[ir]# to get the proper value in length units
                    d_φ_2 = (axφ[iφ]-axφ[iφ-1]) * axr[ir]
                    eφ = ( Δp_φ_1/d_φ_1 + Δp_φ_2/d_φ_2) / 2
                else
                    Δp_φ_1 = p[ir ,iφ+1, iz]-p[ir ,iφ, iz]
                    Δp_φ_2 = p[ir ,iφ, iz]-p[ir ,iφ-1, iz]
                    d_φ_1 = (axφ[iφ+1]-axφ[iφ]) * axr[ir]# to get the proper value in length units
                    d_φ_2 = (axφ[iφ]-axφ[iφ-1]) * axr[ir]
                    eφ = ( Δp_φ_1/d_φ_1 + Δp_φ_2/d_φ_2) / 2
                end
                isinf(eφ) || isnan(eφ) ? eφ = 0.0 : nothing # for small radii and small distances(center of the grid) it would yield Infs or Nans
                if iz-1<1
                    Δp_z_1 = p[ir ,iφ, iz+1]-p[ir ,iφ, iz]
                    d_z_1 = axz[iz+1]-axz[iz]
                    ez = ( Δp_z_1/d_z_1 )
                elseif  iz+1 > size(ef,3)
                    Δp_z_1 = p[ir ,iφ, iz]-p[ir ,iφ, iz-1]
                    d_z_1 = axz[iz]-axz[iz-1]
                    ez = ( Δp_z_1/d_z_1 )
                else
                    Δp_z_1 = p[ir ,iφ, iz+1]-p[ir ,iφ, iz]
                    Δp_z_2 = p[ir ,iφ, iz]-p[ir ,iφ, iz-1]
                    d_z_1 = axz[iz+1]-axz[iz]
                    d_z_2 = axz[iz]-axz[iz-1]
                    ez = ( Δp_z_1/d_z_1 + Δp_z_2/d_z_2) / 2
                end
                if pointtypes[ir, iφ, iz] & update_bit == 0 # boundary points
                    if (1 < ir < size(pointtypes, 1))
                        if (pointtypes[ir - 1, iφ, iz] & update_bit > 0) && (pointtypes[ir + 1, iφ, iz] & update_bit > 0)
                            er = 0
                        elseif (pointtypes[ir - 1, iφ, iz] & update_bit > 0) || (pointtypes[ir + 1, iφ, iz] & update_bit > 0)
                            er *= 2
                        end
                    end
                    if (1 < iφ < size(pointtypes, 2))
                        if (pointtypes[ir, iφ - 1, iz] & update_bit > 0) && (pointtypes[ir, iφ + 1, iz] & update_bit > 0)
                            eφ = 0
                        elseif (pointtypes[ir, iφ - 1, iz] & update_bit > 0) || (pointtypes[ir, iφ + 1, iz] & update_bit > 0)
                            eφ *= 2
                        end
                    end
                    if (1 < iz < size(pointtypes, 3))
                        if (pointtypes[ir, iφ, iz - 1] & update_bit > 0) && (pointtypes[ir, iφ, iz + 1] & update_bit > 0)
                            ez = 0
                        elseif (pointtypes[ir, iφ, iz - 1] & update_bit > 0) || (pointtypes[ir, iφ, iz + 1] & update_bit > 0)
                            ez *= 2
                        end
                    end

                end
                ef[ir,iφ,iz] = [-er, -eφ, -ez]
            end
        end
    end
    if fieldvector_coordinates == :xyz
        ef = convert_field_vectors_to_xyz(ef, axφ)
    end
    return ef
end

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
function get_xyz_vector_from_rφz_vector(v::AbstractArray)::AbstractArray
    return [v[1]*cos(v[2]),v[1]*sin(v[2]),v[3]]
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
#
# function get_xyz_vector_from_field_vector(vector, α)
#         Rα = @SArray([cos(α) -1*sin(α) 0;sin(α) cos(α) 0;0 0 1])
#         result = Rα*vector
#         result
# end
function interpolated_scalarfield(ep::ScalarPotential{T}) where {T}
    knots = ep.grid.axes#(grid.r, grid.φ, grid.z)
    i = interpolate(knots, ep.data, Gridded(Linear()))
    vector_field_itp = extrapolate(i, Periodic())
    return vector_field_itp
end

function interpolated_vectorfield(vectorfield::AbstractArray{<:SVector{3, T},3}, ep::ElectricPotential{T}) where {T}
    knots = ep.grid.axes#(grid.r, grid.φ, grid.z)
    i = interpolate(knots, vectorfield, Gridded(Linear()))
    vectorfield_itp = extrapolate(i, Periodic())
    return vectorfield_itp
end

function setup_interpolated_vectorfield(vectorfield, grid::CylindricalGrid{T}) where {T}
    knots = grid.axes #(grid.r, grid.φ, grid.z)
    i = interpolate(knots, vectorfield, Gridded(Linear()))
    vector_field_itp = extrapolate(i, Periodic())
    return vector_field_itp
end

function get_interpolated_drift_field(velocity_field, grid::CylindricalGrid{T}) where {T}
    knots = grid.axes
    i = interpolate(knots, velocity_field, Gridded(Linear()))
    velocity_field_itp = extrapolate(i, Periodic())
    return velocity_field_itp
end

include("plot_recipes.jl")

function get_electric_field_from_potential(ep::ElectricPotential{T, 3, :Cartesian}, pointtypes::PointTypes{T})::Array{SArray{Tuple{3},T,1,3}, 3} where {T <: AbstractFloat}
    error("Not yet implemented.")
    
    axr::Vector{T} = collect(ep.grid[:x])
    axφ::Vector{T} = collect(ep.grid[:y])
    axz::Vector{T} = collect(ep.grid[:z])
    axr_ext::Vector{T} = get_extended_ticks(ep.grid[:x])
    axφ_ext::Vector{T} = get_extended_ticks(ep.grid[:y])
    axz_ext::Vector{T} = get_extended_ticks(ep.grid[:z])

    ef::Array{SVector{3, T}} = Array{SVector{3, T}}(undef, size(ep.data))


    for iz in eachindex(axz)
        for iz in eachindex(axz)
            for iz in eachindex(axz)

            end
        end
    end
    return ef
end

