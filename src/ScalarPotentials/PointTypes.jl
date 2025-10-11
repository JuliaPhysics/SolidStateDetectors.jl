"""
    const PointType = UInt8

Stores certain information about a grid point via bit-flags. 

Right now, there are:

    const update_bit           = 0x01
    const undepleted_bit       = 0x02
    const pn_junction_bit      = 0x04
    const bulk_bit             = 0x08
    const inactive_layer_bit   = 0x10
    const inactive_contact_bit = 0x20

## Examples

How to get information out of a `PointType` variable `point_type`:
1. `point_type & update_bit == 0` -> do not update this point (for fixed points)     
2. `point_type & update_bit >  0` -> do update this point    
3. `point_type & undepleted_bit > 0` -> this point is undepleted
4. `point_type & pn_junction_bit > 0` -> this point belongs to the solid state detector, meaning that it is in the volume of the n-type or p-type material.
5. `point_type & bulk_bit > 0` -> this point is only surrounded by points marked as `pn_junction_bit`
6. `point_type & inactive_layer_bit > 0` -> this point is part of the inactive layer
7. `point_type & inactive_contact_bit > 0` -> this point is part of the contact next to the inactive layer
"""
const PointType       = UInt8

# Point types for electric potential calculation
const update_bit           = 0x01 # parse(UInt8, "00000001", base=2) # 1 -> do update; 0 -> do not update
const undepleted_bit       = 0x02 # parse(UInt8, "00000010", base=2) # 0 -> depleted point; 1 -> undepleted point
const pn_junction_bit      = 0x04 # parse(UInt8, "00000100", base=2) # 0 -> point is not part of pn-junction; 1 -> point is part of the pn-junction
const bulk_bit             = 0x08 # parse(UInt8, "00001000", base=2) # 0 -> point is surrounded by points that do not belong to the pn-junction; 1 -> point is only surrounded by points in the pn-junction
const inactive_layer_bit   = 0x10 # parse(UInt8, "00010000", base=2) # 0 -> point is not part of the inactive layer; 1 -> point is part of the inactive layer
const inactive_contact_bit = 0x20 # parse(UInt8, "00100000", base=2) # 0 -> point is not part of the contact next to the inactive layer; 1 -> point is part of the contact next to the inactive layer 

is_pn_junction_point_type(p::PointType) = p & pn_junction_bit > 0
is_undepleted_point_type(p::PointType) = p & undepleted_bit > 0
is_fixed_point_type(p::PointType) = p & update_bit == 0
is_in_inactive_layer(p::PointType) = p & inactive_layer_bit > 0
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
@inline getindex(point_types::PointTypes{T, N, S}, pt::AbstractCoordinatePoint{T}) where {T, N, S} = point_types.data[find_closest_gridpoint(pt, point_types.grid)...]

function in(pt::AbstractCoordinatePoint{T}, point_types::PointTypes{T, 3, S})::Bool where {T <: SSDFloat, S}
    cpt = _convert_point(pt, S)
    point_type::PointType = point_types[cpt]
    return point_type & bulk_bit > 0
end


function mark_bulk_bits!(point_types::Array{PointType, 3})
    i1max, i2max, i3max = size(point_types)
    for i1 in Base.OneTo(i1max), i2 in Base.OneTo(i2max), i3 in Base.OneTo(i3max)
        (point_types[i1,i2,i3] & pn_junction_bit == 0) && continue 
        point_types[i1,i2,i3] += bulk_bit * begin
            in_bulk = true
            for j1 in max(i1-1,1):min(i1+1,i1max), j2 in max(i2-1,1):min(i2+1,i2max), j3 in max(i3-1,1):min(i3+1,i3max)
                in_bulk &= (point_types[j1,j2,j3] & pn_junction_bit > 0) &&  (point_types[j1,j2,j3] & update_bit > 0)
                !in_bulk && break
            end
            in_bulk
        end
    end
    point_types
end


"""
    is_depleted(point_types::PointTypes)::Bool
    
Returns `true` if all [`PointType`](@ref) values of
the [`PointTypes`](@ref) of a [`Simulation`](@ref) are marked as depleted
and `false` if any point in the [`PointTypes`](@ref) is marked as undepleted
excluding the inactive layer.

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
    !any(b -> (bulk_bit & b > 0) && (undepleted_bit & b > 0) &&
    (inactive_layer_bit & b == 0), point_types.data)

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
                        active_volume += if iφ == length(point_types.grid.φ) || iφ == 1
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

function Base.NamedTuple(point_types::PointTypes{T, 3}) where {T}
    return (
        grid = NamedTuple(point_types.grid),
        values = point_types.data,
    )
end

Base.convert(T::Type{NamedTuple}, x::PointTypes) = T(x)

"""
    get_inactivelayer_indices(vec::AbstractVector{PointType})::Vector{Vector{Int}}

Scans a 1D vector and returns index groups adjacent to points flagged with `inactive_contact_bit` (points that belongs yo the contact next to the inactive layer). Each group consists of consecutive indices where all of the `required_bits` (`update_bit | undepleted_bit | pn_junction_bit`) are set (undepleted points that belong to the inactive layer).

Note: For detectors with undepleted regions touching the inactive layer, this method is not recommended since it won't identify the inactive layer. In this extreme case, define the inactive layer geometry boundaries via inactive_layer_geometry.
"""

function get_inactivelayer_indices(vec::AbstractVector{PointType})
    pt_indices = Vector{Vector{Int}}()
    len = length(vec)
    required_bits = update_bit | undepleted_bit | pn_junction_bit

    i = 1
    while i <= len
        if (vec[i] & inactive_contact_bit) != 0
            # start checking to the right
            region = Int[]
            j = i + 1
            while j <= len && (vec[j] & required_bits) == required_bits
                push!(region, j)
                j += 1
            end

            # start checking to the left
            j = i - 1
            while j >= 1 && (vec[j] & required_bits) == required_bits
                push!(region, j)
                j -= 1
            end

            if !isempty(region)
                push!(pt_indices, sort(region))
            end
        end
        i += 1
    end

    return pt_indices
end

"""
    mark_inactivelayer_bits!(point_types::Array{PointType, 3})

Applies the bitwise marking rule to slices across all three dimensions of
a 3D array, modifying it in-place.
"""
function mark_inactivelayer_bits!(point_types::Array{PointType, 3})
    sz1, sz2, sz3 = size(point_types)

    # dim 1 (vary i, fix j & k)
    for j in 1:sz2, k in 1:sz3
        vec = @view point_types[:, j, k]
        matches = get_inactivelayer_indices(vec)
        for match in matches, idx in match
            vec[idx] |= inactive_layer_bit
        end
    end

    # dim 2 (vary j, fix i & k)
    for i in 1:sz1, k in 1:sz3
        vec = @view point_types[i, :, k]
        matches = get_inactivelayer_indices(vec)
        for match in matches, idx in match
            vec[idx] |= inactive_layer_bit
        end
    end

    # dim 3 (vary k, fix i & j)
    for i in 1:sz1, j in 1:sz2
        vec = @view point_types[i, j, :]
        matches = get_inactivelayer_indices(vec)
        for match in matches, idx in match
            vec[idx] |= inactive_layer_bit
        end
    end
end


"""
        in_inactive_layer(pt::CartesianPoint{T},
           g::AbstractGeometry{T},
           point_types::PointTypes{T})

Returns if a CartesianPoint belongs to the inactive layer using the geometry of the inactive layer

        in_inactive_layer(pt::CartesianPoint{T},
           ::Nothing,
           point_types::PointTypes{T})
Returns if a CartesianPoint belongs to the inactive layer using point types
"""

@inline in_inactive_layer(pt::CartesianPoint{T}, g::AbstractGeometry{T}, ::PointTypes{T}) where {T} = in(pt, g)
@inline in_inactive_layer(pt::CartesianPoint{T}, ::Nothing, point_types::PointTypes{T}) where {T} = is_in_inactive_layer(point_types[pt])
