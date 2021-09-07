"""
    const PointType = UInt8

Stores certain information about a grid point via bit-flags. 

Right now, there are:

    const update_bit      = 0x01
    const undepleted_bit  = 0x02
    const pn_junction_bit = 0x04

## Examples

How to get information out of a `PointType` variable `point_type`:
1. `point_type & update_bit == 0` -> do not update this point (for fixed points)     
2. `point_type & update_bit >  0` -> do update this point    
3. `point_type & undepleted_bit > 0` -> this point is undepleted
4. `point_type & pn_junction_bit > 0` -> this point belongs to the solid state detector. So it is in the volume of the n-type or p-type material.
"""
const PointType       = UInt8

# Point types for electric potential calculation
const update_bit      = 0x01 # parse(UInt8, "00000001", base=2) # 1 -> do update; 0 -> do not update
const undepleted_bit  = 0x02 # parse(UInt8, "00000010", base=2) # 0 -> depleted point; 1 -> undepleted point
const pn_junction_bit = 0x04 # parse(UInt8, "00001000", base=2) # 0 -> point is not inside a bubble; 1 -> point is inside a bubble

is_pn_junction_point_type(p::PointType) = p & pn_junction_bit > 0
is_undepleted_point_type(p::PointType) = p & undepleted_bit > 0
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
* `grid::Grid{T, N, S, AT}`: [`Grid`](@ref) defining the discrete points for which the point types are determined.

See also [`PointType`](@ref).
"""
struct PointTypes{T, N, S, AT} <: AbstractArray{T, N}
    data::Array{PointType, N}
    grid::Grid{T, N, S, AT}
end

@inline size(point_types::PointTypes{T, N, S}) where {T, N, S} = size(point_types.data)
@inline length(point_types::PointTypes{T, N, S}) where {T, N, S} = length(point_types.data)
@inline getindex(point_types::PointTypes{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(point_types.data, I...)
@inline getindex(point_types::PointTypes{T, N, S}, i::Int) where {T, N, S} = getindex(point_types.data, i)
@inline getindex(point_types::PointTypes{T, N, S}, s::Symbol) where {T, N, S} = getindex(point_types.grid, s)

function in(pt::AbstractCoordinatePoint{T}, point_types::PointTypes{T, 3, S})::Bool where {T <: SSDFloat, S}
    cpt = _convert_point(pt, S)
    grid::Grid{T, 3, S} = point_types.grid
    i1::Int = searchsortednearest(grid.axes[1].ticks, cpt[1])
    i2::Int = searchsortednearest(grid.axes[2].ticks, cpt[2])
    i3::Int = searchsortednearest(grid.axes[3].ticks, cpt[3])
    return point_types.data[i1, i2, i3] & pn_junction_bit > 0
end

"""
    is_depleted(point_types::PointTypes)::Bool
    
Returns `true` if all [`PointType`](@ref) values of
the [`PointTypes`](@ref) of a [`Simulation`](@ref) are marked as depleted
and `false` if any point in the [`PointTypes`](@ref) is marked as undepleted.

It can be used to determine whether the [`SolidStateDetector`](@ref) is
depleted at the provided bias voltage.

## Arguments
* `point_types::PointTypes`: [`PointTypes`](@ref) of a [`Simulation`](@ref).

## Examples
```julia
is_depleted(sim.point_types)
```
"""
is_depleted(point_types::PointTypes)::Bool = 
    !any(b -> undepleted_bit & b > 0, point_types.data)


"""
    get_active_volume(point_types::PointTypes{T}) where {T}

Returns an approximation of the active volume of the detector by summing up the cell volumes of
all cells marked as depleted.

## Arguments
* `point_types::PointTypes{T}`: Point types of a [`Simulation`](@ref).

## Examples
    get_active_volume(sim.point_types)

!!! note
    Only `φ`-symmetries are taken into account.
"""
function get_active_volume(point_types::PointTypes{T, 3, Cylindrical}) where {T}
    active_volume::T = 0

    only_2d::Bool = length(point_types.grid.axes[2]) == 1
    cyclic::T = point_types.grid.axes[2].interval.right

    r_ext::Vector{T} = get_extended_ticks(point_types.grid.axes[1])
    φ_ext::Vector{T} = get_extended_ticks(point_types.grid.axes[2])
    z_ext::Vector{T} = get_extended_ticks(point_types.grid.axes[3])

    mpr::Vector{T} = midpoints(r_ext)
    mpφ::Vector{T} = midpoints(φ_ext)
    mpz::Vector{T} = midpoints(z_ext)
    Δmpφ::Vector{T} = diff(mpφ)
    Δmpz::Vector{T} = diff(mpz)
    Δmpr_squared::Vector{T} = T(0.5) .* ((mpr[2:end].^2) .- (mpr[1:end-1].^2))
    if r_ext[2] == 0
        Δmpr_squared[1] = T(0.5) * (mpr[2]^2)
    end

    isclosed::Bool = typeof(point_types.grid.axes[2].interval).parameters[2] == :closed 
    for iz in eachindex(point_types.grid.axes[3])
        if !isclosed || only_2d
            for iφ in eachindex(point_types.grid.axes[2])
                for ir in eachindex(point_types.grid.axes[1])
                    point_type::PointType = point_types[ir, iφ, iz]
                    if (point_type & pn_junction_bit > 0) && (point_type & undepleted_bit == 0) && (point_type & update_bit > 0)
                        dV::T = Δmpz[iz] * Δmpφ[iφ] * Δmpr_squared[ir]
                        active_volume += dV
                    end
                end
            end
        elseif isclosed && !only_2d
            for iφ in eachindex(point_types.grid.axes[2])
                for ir in eachindex(point_types.grid.axes[1])
                    point_type::PointType = point_types[ir, iφ, iz]
                    if (point_type & pn_junction_bit > 0) && (point_type & undepleted_bit == 0) && (point_type & update_bit > 0)
                        dV = Δmpz[iz] * Δmpφ[iφ] * Δmpr_squared[ir]
                        active_volume += if iφ == length(point_types.φ) || iφ == 1
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

function get_active_volume(point_types::PointTypes{T, 3, Cartesian}) where {T}
    active_volume::T = 0

    x_ext::Vector{T} = get_extended_ticks(point_types.grid.axes[1])
    y_ext::Vector{T} = get_extended_ticks(point_types.grid.axes[2])
    z_ext::Vector{T} = get_extended_ticks(point_types.grid.axes[3])

    mpx::Vector{T} = midpoints(x_ext)
    mpy::Vector{T} = midpoints(y_ext)
    mpz::Vector{T} = midpoints(z_ext)
    Δmpx::Vector{T} = diff(mpx)
    Δmpy::Vector{T} = diff(mpy)
    Δmpz::Vector{T} = diff(mpz)

    for iz in eachindex(point_types.grid.axes[3])
        for iy in eachindex(point_types.grid.axes[2])
            for ix in eachindex(point_types.grid.axes[1])
                point_type::PointType = point_types[ix, iy, iz]
                if (point_type & pn_junction_bit > 0) && (point_type & undepleted_bit == 0) && (point_type & update_bit > 0)
                    dV = Δmpx[ix] * Δmpy[iy] * Δmpz[iz]
                    active_volume += dV
                end
            end
        end
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

function NamedTuple(point_types::PointTypes{T, 3}) where {T}
    return (
        grid = NamedTuple(point_types.grid),
        values = point_types.data,
    )
end

Base.convert(T::Type{NamedTuple}, x::PointTypes) = T(x)


