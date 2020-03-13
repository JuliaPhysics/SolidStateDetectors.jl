"""
    const PointType = UInt8

Stores certain information about a grid point via bit-flags. 

Right now there are:

    `const update_bit      = 0x01`
    `const undepleted_bit  = 0x02`
    `const pn_junction_bit = 0x04`

How to get information out of a PointType variable `pt`:
1. `pt & update_bit == 0` -> do not update this point (for fixed points)     
2. `pt & update_bit >  0` -> do update this point    
3. `pt & undepleted_bit > 0` -> this point is undepleted
4. `pt & pn_junction_bit > 0` -> this point belongs to the solid state detector. So it is in the volume of the n-type or p-type material.
"""
const PointType       = UInt8

# Point types for electric potential calculation
const update_bit      = 0x01 # parse(UInt8, "00000001", base=2) # 1 -> do update; 0 -> do not update
const undepleted_bit  = 0x02 # parse(UInt8, "00000010", base=2) # 0 -> depleted point; 1 -> undepleteded point
const pn_junction_bit = 0x04 # parse(UInt8, "00001000", base=2) # 0 -> point is not inside a bubble; 1 -> point is inside a bubble
# const bubble_bit      = 0x08 # parse(UInt8, "00000100", base=2) # 0 -> point is not inside a bubble; 1 -> point is inside a bubble

const max_pointtype_value = update_bit + undepleted_bit + pn_junction_bit #+ bubble_bit

"""
    PointTypes{T, N, S} <: AbstractArray{T, N}

PointTypes stores the point type of each grid point.
"""
struct PointTypes{T, N, S} <: AbstractArray{T, N}
    data::Array{PointType, N}
    grid::Grid{T, N, S}
end

@inline size(pts::PointTypes{T, N, S}) where {T, N, S} = size(pts.data)
@inline length(pts::PointTypes{T, N, S}) where {T, N, S} = length(pts.data)
@inline getindex(pts::PointTypes{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(pts.data, I...)
@inline getindex(pts::PointTypes{T, N, S}, i::Int) where {T, N, S} = getindex(pts.data, i)
@inline getindex(pts::PointTypes{T, N, S}, s::Symbol) where {T, N, S} = getindex(pts.grid, s)

function in(pt::AbstractCoordinatePoint{T}, pts::PointTypes{T, 3, S})::Bool where {T <: SSDFloat, S}
    i1::Int, i2::Int, i3::Int = searchsortednearest(pt, pts.grid)
    return pts.data[i1, i2, i3] & pn_junction_bit > 0
end

is_depleted(point_types::PointTypes)::Bool = 
    !any(b -> undepleted_bit & b > 0, point_types.data)


"""
    get_active_volume(pts::PointTypes{T}) where {T}
Returns an approximation of the active volume of the detector by summing up the cell volumes of
all depleted cells.
"""
function get_active_volume(pts::PointTypes{T, 3, :cylindrical}) where {T}
    active_volume::T = 0

    only_2d::Bool = length(pts.grid.φ) == 1
    cyclic::T = pts.grid.φ.interval.right

    r_ext::Vector{T} = get_extended_ticks(pts.grid.r)
    φ_ext::Vector{T} = get_extended_ticks(pts.grid.φ)
    z_ext::Vector{T} = get_extended_ticks(pts.grid.z)
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

    isclosed::Bool = typeof(pts.grid.φ.interval).parameters[2] == :closed 
    for iz in eachindex(pts.grid.z)
        if !isclosed || only_2d
            for iφ in eachindex(pts.grid.φ)
                for ir in eachindex(pts.grid.r)
                    pt::PointType = pts[ir, iφ, iz]
                    if (pt & pn_junction_bit > 0) && (pt & undepleted_bit == 0) && (pt & update_bit > 0)
                        dV::T = Δmpz[iz] * Δmpφ[iφ] * Δmpr_squared[ir]
                        active_volume += dV
                    end
                end
            end
        elseif isclosed && !only_2d
            for iφ in eachindex(pts.grid.φ)
                for ir in eachindex(pts.grid.r)
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
    T = typeof(grid[1].ticks[1])
    S = get_coordinate_system(grid)
    N = get_number_of_dimensions(grid)
    PointTypes{T, N, S}( nt.values, grid )
end

Base.convert(T::Type{PointTypes}, x::NamedTuple) = T(x)

function NamedTuple(pts::PointTypes{T, 3}) where {T}
    return (
        grid = NamedTuple(pts.grid),
        values = pts.data,
    )
end

Base.convert(T::Type{NamedTuple}, x::PointTypes) = T(x)




@recipe function f( pts::PointTypes{T, 3, :cylindrical};
                    r = missing,
                    φ = missing,
                    z = missing ) where {T}
    g::Grid{T, 3, :cylindrical} = pts.grid
   
    seriescolor --> :viridis
    st --> :heatmap
    aspect_ratio --> 1
    
    cross_section::Symbol, idx::Int = if ismissing(φ) && ismissing(r) && ismissing(z)
        :φ, 1
    elseif !ismissing(φ) && ismissing(r) && ismissing(z)
        φ_rad::T = T(deg2rad(φ))
        while !(g.φ.interval.left <= φ_rad <= g.φ.interval.right)
            if φ_rad > g.φ.interval.right
                φ_rad -= g.φ.interval.right - g.φ.interval.left
            elseif φ_rad < g.φ.interval.left
                φ_rad += g.φ.interval.right - g.φ.interval.left
            end
        end
        :φ, searchsortednearest(g.φ, φ_rad)
    elseif ismissing(φ) && !ismissing(r) && ismissing(z)
        :r, searchsortednearest(g.r, T(r))
    elseif ismissing(φ) && ismissing(r) && !ismissing(z)
        :z, searchsortednearest(g.z, T(z))
    else
        error(ArgumentError, ": Only one of the keywords `r, φ, z` is allowed.")
    end
    value::T = if cross_section == :φ
        g.φ[idx]
    elseif cross_section == :r    
        g.r[idx]
    elseif cross_section == :z
        g.z[idx]
    end
    
    @series begin
        clims --> (0, max_pointtype_value)
        if cross_section == :φ
            title --> "Point Type Map @$(cross_section) = $(round(rad2deg(value), sigdigits = 2))"
            xlabel --> "r / m"
            ylabel --> "z / m"
            size --> ( 400, 350 / (g.r[end] - g.r[1]) * (g.z[end] - g.z[1]) )
            g.r, g.z, pts.data[:, idx,:]'
        elseif cross_section == :r
            title --> "Point Type Map @$(cross_section) = $(round(value, sigdigits = 2))"
            g.φ, g.z, pts.data[idx,:,:]'
        elseif cross_section == :z
            title --> "Point Type Map @$(cross_section) = $(round(value, sigdigits = 2))"
            proj --> :polar
            g.φ, g.r, pts.data[:,:,idx]
        end
    end
end



@recipe function f( pts::PointTypes{T, 3, :cartesian};
                    x = missing,
                    y = missing,
                    z = missing ) where {T}
    g::Grid{T, 3, :cartesian} = pts.grid
   
    seriescolor --> :viridis
    st --> :heatmap
    aspect_ratio --> 1
    
    cross_section::Symbol, idx::Int = if ismissing(x) && ismissing(y) && ismissing(z)
        :x, 1
    elseif !ismissing(x) && ismissing(y) && ismissing(z)
        :x, searchsortednearest(g[:x], T(x))
    elseif ismissing(x) && !ismissing(y) && ismissing(z)
        :y, searchsortednearest(g[:y], T(y))
    elseif ismissing(x) && ismissing(y) && !ismissing(z)
        :z, searchsortednearest(g.z, T(z))
    else
        error(ArgumentError, ": Only one of the keywords `r, y, z` is allowed.")
    end
    value::T = if cross_section == :x
        g[:x][idx]
    elseif cross_section == :y
        g[:y][idx]
    elseif cross_section == :z
        g.z[idx]
    end
    
    @series begin
        clims --> (0, max_pointtype_value)
        title --> "Point Type Map @$(cross_section) = $(round(value, sigdigits = 2))"
        if cross_section == :x
            xlabel --> "y / m"
            ylabel --> "z / m"
            g[:y], g.z, pts.data[idx, :, :]'
        elseif cross_section == :y
            xlabel --> "x / m"
            ylabel --> "z / m"
            g[:x], g.z, pts.data[:, idx, :]'
        elseif cross_section == :z
            xlabel --> "x / m"
            ylabel --> "y / m"
            g[:x], g[:y], pts.data[:,:,idx]'
        end
    end
end
