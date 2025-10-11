"""
    struct ElectricField{T, N, S, AT} <: AbstractArray{T, N}
        
Electric field of the simulation in units of volt per meter (V/m).
        
## Parametric types 
* `T`: Element type of `grid.axes`.
* `N`: Dimension of the `grid` and `data` array.  
* `S`: Coordinate system (`Cartesian` or `Cylindrical`).
* `AT`: Axes type.
        
## Fields
* `data::Array{<:StaticArray{Tuple{N}, T}, N}`: Array containing the field vectors of the electric field at the discrete points of the `grid`.
* `grid::Grid{T, N, S, AT}`: [`Grid`](@ref) defining the discrete points for which the electric field is determined.
"""
struct ElectricField{T, N, S, AT} <: AbstractArray{T, N}
    data::Array{<:StaticArray{Tuple{N}, T}, N}
    grid::Grid{T, N, S, AT}
end

@inline size(ef::ElectricField{T, N, S}) where {T, N, S} = size(ef.data)
@inline length(ef::ElectricField{T, N, S}) where {T, N, S} = length(ef.data)
@inline getindex(ef::ElectricField{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(ef.data, I...)
@inline getindex(ef::ElectricField{T, N, S}, i::Int) where {T, N, S} = getindex(ef.data, i)
@inline getindex(ef::ElectricField{T, N, S}, s::Symbol) where {T, N, S} = getindex(ef.grid, s)


function Base.NamedTuple(ef::ElectricField{T, 3}) where {T}
    return (
        grid = NamedTuple(ef.grid),
        values = ef.data #* internal_efield_unit,
    )
end
Base.convert(T::Type{NamedTuple}, x::ElectricField) = T(x)

function ElectricField(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = eltype(grid)
    S = get_coordinate_system(grid)
    N = get_number_of_dimensions(grid)
    ef = Array{SVector{3, T}}(undef, size(grid)...)
    @inbounds for i in eachindex(ef)
        # ef[i] = ustrip.(uconvert.(internal_efield_unit, nt.values[i]))
        ef[i] = nt.values[i]
    end
    ElectricField{T, N, S, typeof(grid.axes)}( ef, grid)
end
Base.convert(T::Type{ElectricField}, x::NamedTuple) = T(x)



function ElectricField(epot::ElectricPotential{T, 3, S}, point_types::PointTypes{T}; use_nthreads::Int = Base.Threads.nthreads()) where {T, S}
    return ElectricField{T, 3, S, typeof(grid.axes)}(get_electric_field_from_potential( epot, point_types; use_nthreads ), epot.grid)
end


function get_electric_field_from_potential(epot::ElectricPotential{T, 3, Cylindrical}, point_types::PointTypes{T}; use_nthreads::Int = Base.Threads.threadid())::ElectricField{T, 3, Cylindrical} where {T <: SSDFloat}
    p = epot.data
    axr::Vector{T} = collect(epot.grid.axes[1])
    axφ::Vector{T} = collect(epot.grid.axes[2])
    axz::Vector{T} = collect(epot.grid.axes[3])

    cyclic::T = epot.grid.axes[2].interval.right
    ef = Array{SVector{3, T}}(undef, size(p)...)
    @onthreads 1:use_nthreads for iz in workpart(1:size(ef, 3), 1:use_nthreads, Base.Threads.threadid())
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
                    d_φ_2 = (cyclic + axφ[iφ] - axφ[end]) * axr[ir]
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
                if point_types[ir, iφ, iz] & update_bit == 0 # boundary points
                    if (1 < ir < size(point_types, 1))
                        if (point_types[ir - 1, iφ, iz] & update_bit > 0) && (point_types[ir + 1, iφ, iz] & update_bit > 0)
                            if (point_types[ir - 1, iφ, iz] & pn_junction_bit > 0) && (point_types[ir + 1, iφ, iz] & pn_junction_bit == 0)
                                er = Δp_r_2 / d_r_2
                            elseif (point_types[ir + 1, iφ, iz] & pn_junction_bit > 0) && (point_types[ir - 1, iφ, iz] & pn_junction_bit == 0)
                                er = Δp_r_1 / d_r_1
                            end
                        elseif (point_types[ir - 1, iφ, iz] & update_bit > 0) || (point_types[ir + 1, iφ, iz] & update_bit > 0)
                            er *= 2
                        end
                    end
                    if (1 < iφ < size(point_types, 2))
                        if (point_types[ir, iφ - 1, iz] & update_bit > 0) && (point_types[ir, iφ + 1, iz] & update_bit > 0)
                            if (point_types[ir, iφ - 1, iz] & pn_junction_bit > 0) && (point_types[ir, iφ + 1, iz] & pn_junction_bit == 0)
                                eφ = Δp_φ_2 / d_φ_2
                            elseif (point_types[ir, iφ + 1, iz] & pn_junction_bit > 0) && (point_types[ir, iφ - 1, iz] & pn_junction_bit == 0)
                                eφ = Δp_φ_1 / d_φ_1
                            end
                        elseif (point_types[ir, iφ - 1, iz] & update_bit > 0) || (point_types[ir, iφ + 1, iz] & update_bit > 0)
                            eφ *= 2
                        end
                    end
                    if (1 < iz < size(point_types, 3))
                        if (point_types[ir, iφ, iz - 1] & update_bit > 0) && (point_types[ir, iφ, iz + 1] & update_bit > 0)
                            if (point_types[ir, iφ, iz - 1] & pn_junction_bit > 0) && (point_types[ir, iφ, iz + 1] & pn_junction_bit == 0)
                                ez = Δp_z_2 / d_z_2
                            elseif (point_types[ir, iφ, iz + 1] & pn_junction_bit > 0) && (point_types[ir, iφ, iz - 1] & pn_junction_bit == 0)
                                ez = Δp_z_1 / d_z_1
                            end
                        elseif (point_types[ir, iφ, iz - 1] & update_bit > 0) || (point_types[ir, iφ, iz + 1] & update_bit > 0)
                            ez *= 2
                        end
                    end

                end
                ef[ir,iφ,iz] = [-er, -eφ, -ez]
            end
        end
    end
    ef = convert_field_vectors_to_xyz(ef, axφ)
    return ElectricField(ef, point_types.grid)
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



function interpolated_scalarfield(spot::ScalarPotential{T, 3, Cylindrical}) where {T}
    @inbounds knots = spot.grid.axes[1].ticks, cat(spot.grid.axes[2].ticks,T(2π),dims=1), spot.grid.axes[3].ticks
    ext_data = cat(spot.data, spot.data[:,1:1,:], dims=2)
    i = interpolate!(knots, ext_data, Gridded(Linear()))
    vector_field_itp = extrapolate(i, (Interpolations.Line(), Periodic(), Interpolations.Line()))
    return vector_field_itp
end
function interpolated_scalarfield(spot::ScalarPotential{T, 3, Cartesian}) where {T}
    @inbounds knots = spot.grid.axes[1].ticks, spot.grid.axes[2].ticks, spot.grid.axes[3].ticks
    i = interpolate!(knots, spot.data, Gridded(Linear()))
    vector_field_itp = extrapolate(i, (Interpolations.Line(), Interpolations.Line(), Interpolations.Line()))
    return vector_field_itp
end


function interpolated_vectorfield(vectorfield, grid::CylindricalGrid{T}) where {T}
    extended_vectorfield = cat(vectorfield, vectorfield[:,1:1,:], dims=2)
    @inbounds knots = grid.axes[1].ticks, cat(grid.axes[2].ticks,T(2π),dims=1), grid.axes[3].ticks
    i = interpolate!(knots, extended_vectorfield, Gridded(Linear()))
    velocity_field_itp = extrapolate(i, (Interpolations.Line(), Periodic(), Interpolations.Line()))
    return velocity_field_itp
end
function interpolated_vectorfield(vectorfield, grid::CartesianGrid3D{T}) where {T}
    @inbounds knots = grid.axes[1].ticks, grid.axes[2].ticks, grid.axes[3].ticks
    i = interpolate!(knots, vectorfield, Gridded(Linear()))
    velocity_field_itp = extrapolate(i, (Interpolations.Line(), Interpolations.Line(), Interpolations.Line()))
    return velocity_field_itp
end

function interpolated_vectorfield(ef::ElectricField)
    interpolated_vectorfield(ef.data, ef.grid)
end



function get_electric_field_from_potential(epot::ElectricPotential{T, 3, Cartesian}, point_types::PointTypes{T}; use_nthreads::Int = Base.Threads.nthreads())::ElectricField{T, 3, Cartesian} where {T <: SSDFloat}
    axx::Vector{T} = collect(epot.grid.axes[1])
    axy::Vector{T} = collect(epot.grid.axes[2])
    axz::Vector{T} = collect(epot.grid.axes[3])
    # axx_ext::Vector{T} = get_extended_ticks(epot.grid.axes[1])
    # axy_ext::Vector{T} = get_extended_ticks(epot.grid.axes[2])
    # axz_ext::Vector{T} = get_extended_ticks(epot.grid.axes[3])

    ef::Array{SVector{3, T}} = Array{SVector{3, T}}(undef, size(epot.data))

    @onthreads 1:use_nthreads for ix in workpart(eachindex(axx), 1:use_nthreads, Base.Threads.threadid())
        for iy in eachindex(axy)
            for iz in eachindex(axz)
                if ix - 1 < 1
                    Δp_x_1::T = epot.data[ix + 1, iy, iz] - epot.data[ix, iy, iz]
                    d_x_1::T = axx[ix + 1] - axx[ix]
                    ex::T =  Δp_x_1 / d_x_1
                elseif ix + 1 > size(ef, 1)
                    Δp_x_1 = epot.data[ix, iy, iz] - epot.data[ix - 1, iy, iz]
                    d_x_1 = axx[ix] - axx[ix - 1]
                    ex = Δp_x_1 / d_x_1
                else
                    Δp_x_1 = epot.data[ix + 1, iy, iz] - epot.data[ix ,iy, iz]
                    Δp_x_2::T = epot.data[ix, iy, iz] - epot.data[ix - 1, iy, iz]
                    d_x_1 = axx[ix + 1] - axx[ix]
                    d_x_2::T = axx[ix] - axx[ix - 1]
                    ex = (Δp_x_1 / d_x_1 + Δp_x_2 / d_x_2) / 2
                end

                if iy - 1 < 1
                    Δp_y_1::T = epot.data[ix, iy + 1, iz] - epot.data[ix ,iy, iz]
                    d_y_1::T = axy[iy + 1] - axy[iy]
                    ey::T =  Δp_y_1 / d_y_1
                elseif iy + 1 > size(ef, 2)
                    Δp_y_1 = epot.data[ix, iy, iz] - epot.data[ix, iy - 1, iz]
                    d_y_1 = axy[iy] - axy[iy - 1]
                    ey = Δp_y_1 / d_y_1
                else
                    Δp_y_1 = epot.data[ix, iy + 1, iz] - epot.data[ix ,iy, iz]
                    Δp_y_2::T = epot.data[ix, iy, iz] - epot.data[ix, iy - 1, iz]
                    d_y_1 = axy[iy + 1] - axy[iy]
                    d_y_2::T = axy[iy] - axy[iy - 1]
                    ey = (Δp_y_1 / d_y_1 + Δp_y_2 / d_y_2) / 2
                end

                if iz - 1 < 1
                    Δp_z_1::T = epot.data[ix, iy, iz + 1] - epot.data[ix, iy, iz]
                    d_z_1::T = axz[iz + 1] - axz[iz]
                    ez::T =  Δp_z_1 / d_z_1
                elseif iz + 1 > size(ef, 3)
                    Δp_z_1 = epot.data[ix, iy, iz] - epot.data[ix, iy, iz - 1]
                    d_z_1 = axz[iz] - axz[iz - 1]
                    ez = Δp_z_1 / d_z_1
                else
                    Δp_z_1 = epot.data[ix, iy, iz + 1] - epot.data[ix ,iy, iz]
                    Δp_z_2::T = epot.data[ix, iy, iz] - epot.data[ix, iy, iz - 1]
                    d_z_1 = axz[iz + 1] - axz[iz]
                    d_z_2::T = axz[iz] - axz[iz - 1]
                    ez = (Δp_z_1 / d_z_1 + Δp_z_2 / d_z_2) / 2
                end

                if point_types[ix, iy, iz] & update_bit == 0 # boundary points
                    if (1 < ix < size(point_types, 1))
                        if (point_types[ix - 1, iy, iz] & update_bit > 0) && (point_types[ix + 1, iy, iz] & update_bit > 0)
                            if (point_types[ix - 1, iy, iz] & pn_junction_bit > 0) && (point_types[ix + 1, iy, iz] & pn_junction_bit == 0)
                                ex = Δp_x_2 / d_x_2
                            elseif (point_types[ix + 1, iy, iz] & pn_junction_bit > 0) && (point_types[ix - 1, iy, iz] & pn_junction_bit == 0)
                                ex = Δp_x_1 / d_x_1
                            end
                        elseif (point_types[ix - 1, iy, iz] & update_bit > 0) || (point_types[ix + 1, iy, iz] & update_bit > 0)
                            ex *= 2
                        end
                    end
                    if (1 < iy < size(point_types, 2))
                        if (point_types[ix, iy - 1, iz] & update_bit > 0) && (point_types[ix, iy + 1, iz] & update_bit > 0)
                            if (point_types[ix, iy - 1, iz] & pn_junction_bit > 0) && (point_types[ix, iy + 1, iz] & pn_junction_bit == 0)
                                ey = Δp_y_2 / d_y_2
                            elseif (point_types[ix, iy + 1, iz] & pn_junction_bit > 0) && (point_types[ix, iy - 1, iz] & pn_junction_bit == 0)
                                ey = Δp_y_1 / d_y_1
                            end
                        elseif (point_types[ix, iy - 1, iz] & update_bit > 0) || (point_types[ix, iy + 1, iz] & update_bit > 0)
                            ey *= 2
                        end
                    end
                    if (1 < iz < size(point_types, 3))
                        if (point_types[ix, iy, iz - 1] & update_bit > 0) && (point_types[ix, iy, iz + 1] & update_bit > 0)
                            if (point_types[ix, iy, iz - 1] & pn_junction_bit > 0) && (point_types[ix, iy, iz + 1] & pn_junction_bit == 0)
                                ez = Δp_z_2 / d_z_2
                            elseif (point_types[ix, iy, iz + 1] & pn_junction_bit > 0) && (point_types[ix, iy, iz - 1] & pn_junction_bit == 0)
                                ez = Δp_z_1 / d_z_1
                            end
                        elseif (point_types[ix, iy, iz - 1] & update_bit > 0) || (point_types[ix, iy, iz + 1] & update_bit > 0)
                            ez *= 2
                        end
                    end
                end
                ef[ix, iy, iz] = @SVector [-ex, -ey, -ez]
            end
        end
    end
    return ElectricField(ef, point_types.grid)
end
