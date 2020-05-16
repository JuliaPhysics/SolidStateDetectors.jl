struct ElectricField{T, N, S} <: AbstractArray{T, N}
    data::Array{<:StaticArray{Tuple{N}, T}, N}
    grid::Grid{T, N, S}
end

@inline size(ep::ElectricField{T, N, S}) where {T, N, S} = size(ep.data)
@inline length(ep::ElectricField{T, N, S}) where {T, N, S} = length(ep.data)
@inline getindex(ep::ElectricField{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(ep.data, I...)
@inline getindex(ep::ElectricField{T, N, S}, i::Int) where {T, N, S} = getindex(ep.data, i)
@inline getindex(ep::ElectricField{T, N, S}, s::Symbol) where {T, N, S} = getindex(ep.grid, s)


function NamedTuple(ep::ElectricField{T, 3}) where {T}
    return (
        grid = NamedTuple(ep.grid),
        values = ep.data #* internal_efield_unit,
    )
end
Base.convert(T::Type{NamedTuple}, x::ElectricField) = T(x)

function ElectricField(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = eltype(ustrip.(nt.values[1]))
    S = get_coordinate_system(grid)
    N = get_number_of_dimensions(grid)
    ef = Array{SVector{3, T}}(undef, size(grid)...)
    @inbounds for i in eachindex(ef)
        # ef[i] = ustrip.(uconvert.(internal_efield_unit, nt.values[i]))
        ef[i] = nt.values[i]
    end
    ElectricField{T, N, S}( ef, grid)
end
Base.convert(T::Type{ElectricField}, x::NamedTuple) = T(x)



function ElectricField(ep::ElectricPotential{T, 3, S}, pointtypes::PointTypes{T}) where {T, S}
    return ElectricField{T, 3, S}(get_electric_field_from_potential( ep, pointtypes ), ep.grid)
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


function get_electric_field_from_potential(ep::ElectricPotential{T, 3, :cylindrical}, pointtypes::PointTypes{T}, fieldvector_coordinates=:xyz)::ElectricField{T, 3, :cylindrical} where {T <: SSDFloat}
    p = ep.data
    axr::Vector{T} = collect(ep.grid.axes[1])
    axφ::Vector{T} = collect(ep.grid.axes[2])
    axz::Vector{T} = collect(ep.grid.axes[3])

    cyclic::T = ep.grid.axes[2].interval.right
    ef = Array{SVector{3, T}}(undef, size(p)...)
    for iz in 1:size(ef, 3)
        for iφ in 1:size(ef, 2)
            for ir in 1:size(ef, 1)
                ### r ###
                if ir == 1
                    Δp_r_1 = p[ir+1 ,iφ, iz] - p[ir ,iφ, iz]
                    d_r_1 = axr[ir+1] - axr[ir]
                    er = Δp_r_1 / d_r_1 
                elseif ir == size(ef,1)
                    Δp_r_1 = p[ir ,iφ, iz] - p[ir-1 ,iφ, iz]
                    d_r_1 = axr[ir] - axr[ir-1]
                    er = Δp_r_1 / d_r_1 
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
    return ElectricField(ef, pointtypes.grid)
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
            end
        end
    end
    return field_xyz
end


function interpolated_scalarfield(ep::ScalarPotential{T, 3, :cylindrical}) where {T}
    @inbounds knots = ep.grid.axes[1].ticks, cat(ep.grid.axes[2].ticks,T(2π),dims=1), ep.grid.axes[3].ticks
    ext_data = cat(ep.data, ep.data[:,1:1,:], dims=2)
    i = interpolate(knots, ext_data, Gridded(Linear()))
    vector_field_itp = extrapolate(i, (Interpolations.Line(), Periodic(), Interpolations.Line()))
    return vector_field_itp
end
function interpolated_scalarfield(ep::ScalarPotential{T, 3, :cartesian}) where {T}
    @inbounds knots = ep.grid.axes[1].ticks, ep.grid.axes[2].ticks, ep.grid.axes[3].ticks
    i = interpolate(knots, ep.data, Gridded(Linear()))
    vector_field_itp = extrapolate(i, (Interpolations.Line(), Interpolations.Line(), Interpolations.Line()))
    return vector_field_itp
end


function get_interpolated_drift_field(velocityfield, grid::CylindricalGrid{T}) where {T}
    extended_velocityfield = cat(velocityfield, velocityfield[:,1:1,:], dims=2)
    @inbounds knots = grid.axes[1].ticks, cat(grid.axes[2].ticks,T(2π),dims=1), grid.axes[3].ticks
    i = interpolate(knots, extended_velocityfield, Gridded(Linear()))
    velocity_field_itp = extrapolate(i, (Interpolations.Line(), Periodic(), Interpolations.Line()))
    return velocity_field_itp
end
function get_interpolated_drift_field(velocityfield, grid::CartesianGrid{T}) where {T}
    @inbounds knots = grid.axes[1].ticks, grid.axes[2].ticks, grid.axes[3].ticks
    i = interpolate(knots, velocityfield, Gridded(Linear()))
    velocity_field_itp = extrapolate(i, (Interpolations.Line(), Interpolations.Line(), Interpolations.Line()))
    return velocity_field_itp
end

include("plot_recipes.jl")

function get_electric_field_from_potential(ep::ElectricPotential{T, 3, :cartesian}, pointtypes::PointTypes{T})::ElectricField{T, 3, :cartesian} where {T <: SSDFloat}
    axx::Vector{T} = collect(ep.grid.axes[1])
    axy::Vector{T} = collect(ep.grid.axes[2])
    axz::Vector{T} = collect(ep.grid.axes[3])
    axx_ext::Vector{T} = get_extended_ticks(ep.grid.axes[1])
    axy_ext::Vector{T} = get_extended_ticks(ep.grid.axes[2])
    axz_ext::Vector{T} = get_extended_ticks(ep.grid.axes[3])

    ef::Array{SVector{3, T}} = Array{SVector{3, T}}(undef, size(ep.data))

    for ix in eachindex(axx)
        for iy in eachindex(axy)
            for iz in eachindex(axz)
                if ix - 1 < 1
                    Δp_x_1::T = ep.data[ix + 1, iy, iz] - ep.data[ix, iy, iz]
                    d_x_1::T = axx[ix + 1] - axx[ix]
                    ex::T =  Δp_x_1 / d_x_1
                elseif ix + 1 > size(ef, 1)
                    Δp_x_1 = ep.data[ix, iy, iz] - ep.data[ix - 1, iy, iz]
                    d_x_1 = axx[ix] - axx[ix - 1]
                    ex = Δp_x_1 / d_x_1
                else
                    Δp_x_1 = ep.data[ix + 1, iy, iz] - ep.data[ix ,iy, iz]
                    Δp_x_2::T = ep.data[ix, iy, iz] - ep.data[ix - 1, iy, iz]
                    d_x_1 = axx[ix + 1] - axx[ix]
                    d_x_2::T = axx[ix] - axx[ix - 1]
                    ex = (Δp_x_1 / d_x_1 + Δp_x_2 / d_x_2) / 2
                end

                if iy - 1 < 1
                    Δp_y_1::T = ep.data[ix, iy + 1, iz] - ep.data[ix ,iy, iz]
                    d_y_1::T = axy[iy + 1] - axy[iy]
                    ey::T =  Δp_y_1 / d_y_1
                elseif iy + 1 > size(ef, 2)
                    Δp_y_1 = ep.data[ix, iy, iz] - ep.data[ix, iy - 1, iz]
                    d_y_1 = axy[iy] - axy[iy - 1]
                    ey = Δp_y_1 / d_y_1
                else
                    Δp_y_1 = ep.data[ix, iy + 1, iz] - ep.data[ix ,iy, iz]
                    Δp_y_2::T = ep.data[ix, iy, iz] - ep.data[ix, iy - 1, iz]
                    d_y_1 = axy[iy + 1] - axy[iy]
                    d_y_2::T = axy[iy] - axy[iy - 1]
                    ey = (Δp_y_1 / d_y_1 + Δp_y_2 / d_y_2) / 2
                end

                if iz - 1 < 1
                    Δp_z_1::T = ep.data[ix, iy, iz + 1] - ep.data[ix, iy, iz]
                    d_z_1::T = axz[iz + 1] - axz[iz]
                    ez::T =  Δp_z_1 / d_z_1
                elseif iz + 1 > size(ef, 3)
                    Δp_z_1 = ep.data[ix, iy, iz] - ep.data[ix, iy, iz - 1]
                    d_z_1 = axz[iz] - axz[iz - 1]
                    ez = Δp_z_1 / d_z_1
                else
                    Δp_z_1 = ep.data[ix, iy, iz + 1] - ep.data[ix ,iy, iz]
                    Δp_z_2::T = ep.data[ix, iy, iz] - ep.data[ix, iy, iz - 1]
                    d_z_1 = axz[iz + 1] - axz[iz]
                    d_z_2::T = axz[iz] - axz[iz - 1]
                    ez = (Δp_z_1 / d_z_1 + Δp_z_2 / d_z_2) / 2
                end

                if pointtypes[ix, iy, iz] & update_bit == 0 # boundary points
                    if (1 < ix < size(pointtypes, 1))
                        if (pointtypes[ix - 1, iy, iz] & update_bit > 0) && (pointtypes[ix + 1, iy, iz] & update_bit > 0)
                            ex = 0
                        elseif (pointtypes[ix - 1, iy, iz] & update_bit > 0) || (pointtypes[ix + 1, iy, iz] & update_bit > 0)
                            ex *= 2
                        end
                    end
                    if (1 < iy < size(pointtypes, 2))
                        if (pointtypes[ix, iy - 1, iz] & update_bit > 0) && (pointtypes[ix, iy + 1, iz] & update_bit > 0)
                            ey = 0
                        elseif (pointtypes[ix, iy - 1, iz] & update_bit > 0) || (pointtypes[ix, iy + 1, iz] & update_bit > 0)
                            ey *= 2
                        end
                    end
                    if (1 < iz < size(pointtypes, 3))
                        if (pointtypes[ix, iy, iz - 1] & update_bit > 0) && (pointtypes[ix, iy, iz + 1] & update_bit > 0)
                            ez = 0
                        elseif (pointtypes[ix, iy, iz - 1] & update_bit > 0) || (pointtypes[ix, iy, iz + 1] & update_bit > 0)
                            ez *= 2
                        end
                    end
                end
                ef[ix, iy, iz] = @SVector [-ex, -ey, -ez]
            end
        end
    end
    return ElectricField(ef, pointtypes.grid)
end

function get_electric_field_strength(ef::ElectricField{T}) where {T <: SSDFloat}
    efs::Array{T, 3} = Array{T, 3}(undef, size(ef.data))
    @inbounds for i in eachindex(ef.data)
        efs[i] = norm(ef.data[i])
    end
    return efs
end
