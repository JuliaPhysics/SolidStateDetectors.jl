"""
    const PointType = UInt8

Stores certain information about a grid point via bit-flags. 

Right now there are:

    const update_bit      = 0x01
    const undepleted_bit  = 0x02
    const pn_junction_bit = 0x04

How to get information out of a PointType variable `pt`:
1. `pt & update_bit == 0` -> do not update this point (for fixed points)     
2. `pt & update_bit >  0` -> do update this point    
3. `pt & undepleted_bit > 0` -> this point is undepleted
4. `pt & pn_junction_bit > 0` -> this point belongs to the solid state detector. So it is in the volume of the n-type or p-type material.
"""
const PointType       = UInt8

# Point types for electric potential calculation
const update_bit      = 0x01 # parse(UInt8, "00000001", base=2) # 1 -> do update; 0 -> do not update
const undepleted_bit  = 0x02 # parse(UInt8, "00000010", base=2) # 0 -> depleted point; 1 -> undepleted point
const pn_junction_bit = 0x04 # parse(UInt8, "00001000", base=2) # 0 -> point is not inside a bubble; 1 -> point is inside a bubble

"""
    struct PointTypes{T, N, S, AT} <: AbstractArray{T, N}
        
Information about the grid points used to calculate the [`ElectricPotential`](@ref)
stored via bit-flags. Data is stored as [`PointType`](@ref) which is an `UInt8`.
        
## Parametric types 
* `T`: Element type of `grid.axes`.
* `N`: Dimension of the `grid` and `data` array.  
* `S`: Coordinate system (`Cartesian` or `Cylindrical`).
* `AT`: Axes type.  
        
## Fields
* `data::Array{PointType, N}`: Array containing the point type values at the discrete points of the `grid`.
* `grid::Grid{T, N, S, AT}`: Grid defining the discrete points for which the point types are determined.

See also [`PointType`](@ref).
"""
struct PointTypes{T, N, S, AT} <: AbstractArray{T, N}
    data::Array{PointType, N}
    grid::Grid{T, N, S, AT}
end

@inline size(pts::PointTypes{T, N, S}) where {T, N, S} = size(pts.data)
@inline length(pts::PointTypes{T, N, S}) where {T, N, S} = length(pts.data)
@inline getindex(pts::PointTypes{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(pts.data, I...)
@inline getindex(pts::PointTypes{T, N, S}, i::Int) where {T, N, S} = getindex(pts.data, i)
@inline getindex(pts::PointTypes{T, N, S}, s::Symbol) where {T, N, S} = getindex(pts.grid, s)

function in(pt::AbstractCoordinatePoint{T}, pts::PointTypes{T, 3, S})::Bool where {T <: SSDFloat, S}
    cpt = _convert_point(pt, S)
    g::Grid{T, 3, S} = pts.grid
    i1::Int = searchsortednearest(g.axes[1].ticks, cpt[1])
    i2::Int = searchsortednearest(g.axes[2].ticks, cpt[2])
    i3::Int = searchsortednearest(g.axes[3].ticks, cpt[3])
    return pts.data[i1, i2, i3] & pn_junction_bit > 0
end

"""
    is_depleted(point_types::PointTypes)::Bool
    
Returns a `Bool` value which is `true` if all point types are marked as depleted
and `false` if any point in the point types is marked as undepleted.

## Examples
    is_depleted(sim.point_types)
"""
is_depleted(point_types::PointTypes)::Bool = 
    !any(b -> undepleted_bit & b > 0, point_types.data)


"""
    get_active_volume(pts::PointTypes{T}) where {T}
Returns an approximation of the active volume of the detector by summing up the cell volumes of
all depleted cells.
"""
function get_active_volume(pts::PointTypes{T, 3, Cylindrical}) where {T}
    active_volume::T = 0

    only_2d::Bool = length(pts.grid.axes[2]) == 1
    cyclic::T = pts.grid.axes[2].interval.right

    r_ext::Vector{T} = get_extended_ticks(pts.grid.axes[1])
    φ_ext::Vector{T} = get_extended_ticks(pts.grid.axes[2])
    z_ext::Vector{T} = get_extended_ticks(pts.grid.axes[3])
    Δr_ext::Vector{T} = diff(r_ext)
    Δφ_ext::Vector{T} = diff(φ_ext)
    Δz_ext::Vector{T} = diff(z_ext)

    mpr::Vector{T} = midpoints(r_ext)
    mpφ::Vector{T} = midpoints(φ_ext)
    mpz::Vector{T} = midpoints(z_ext)
    Δmpφ::Vector{T} = diff(mpφ)
    Δmpz::Vector{T} = diff(mpz)
    Δmpr_squared::Vector{T} = T(0.5) .* ((mpr[2:end].^2) .- (mpr[1:end-1].^2))
    if r_ext[2] == 0
        Δmpr_squared[1] = T(0.5) * (mpr[2]^2)
    end

    isclosed::Bool = typeof(pts.grid.axes[2].interval).parameters[2] == :closed 
    for iz in eachindex(pts.grid.axes[3])
        if !isclosed || only_2d
            for iφ in eachindex(pts.grid.axes[2])
                for ir in eachindex(pts.grid.axes[1])
                    pt::PointType = pts[ir, iφ, iz]
                    if (pt & pn_junction_bit > 0) && (pt & undepleted_bit == 0) && (pt & update_bit > 0)
                        dV::T = Δmpz[iz] * Δmpφ[iφ] * Δmpr_squared[ir]
                        active_volume += dV
                    end
                end
            end
        elseif isclosed && !only_2d
            for iφ in eachindex(pts.grid.axes[2])
                for ir in eachindex(pts.grid.axes[1])
                    pt::PointType = pts[ir, iφ, iz]
                    if (pt & pn_junction_bit > 0) && (pt & undepleted_bit == 0) && (pt & update_bit > 0)
                        dV = Δmpz[iz] * Δmpφ[iφ] * Δmpr_squared[ir]
                        active_volume += if iφ == length(pts.φ) || iφ == 1
                            dV / 2
                        else
                            dV
                        end
                    end
                end
            end
        end
    end
    if cyclic > 0
        active_volume *= 2π / cyclic
    end

    f::T = 10^6
    return active_volume * f * Unitful.cm * Unitful.cm * Unitful.cm
end



function PointTypes(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(grid.axes[1].ticks[1])
    S = get_coordinate_system(grid)
    N = get_number_of_dimensions(grid)
    PointTypes{T, N, S, typeof(grid.axes)}( nt.values, grid )
end

Base.convert(T::Type{PointTypes}, x::NamedTuple) = T(x)

function NamedTuple(pts::PointTypes{T, 3}) where {T}
    return (
        grid = NamedTuple(pts.grid),
        values = pts.data,
    )
end

Base.convert(T::Type{NamedTuple}, x::PointTypes) = T(x)


