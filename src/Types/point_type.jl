"""
    const PointType = UInt8

Stores certain information about a grid point via bit-flags. 

Right now there are:

    `const update_bit      = 0x01`
    `const undepleted_bit  = 0x02`
    `const bubble_bit      = 0x04`
    `const pn_junction_bit = 0x08`

How to get information out of a PointType variable `pt`:
1. `pt & update_bit == 0` -> do not update this point (for fixed points)     
2. `pt & update_bit >  0` -> do update this point    
3. `pt & undepleted_bit > 0` -> this point is undepleted
4. `pt & bubble_bit > 0`     -> this point is in an "bubble"-area/volume. So an undepleted region which is not in contact with an electrode (fixed points).
5. `pt & pn_junction_bit > 0`  -> this point to the solid state detector. So it is in the volume of the n-type or p-type material.
"""
const PointType       = UInt8
const update_bit      = 0x01 # parse(UInt8, "00000001", base=2) # 1 -> do update; 0 -> do not update
const undepleted_bit  = 0x02 # parse(UInt8, "00000010", base=2) # 0 -> depleted point; 1 -> undepleteded point
const bubble_bit      = 0x04 # parse(UInt8, "00000100", base=2) # 0 -> point is not inside a bubble; 1 -> point is inside a bubble
const pn_junction_bit = 0x08 # parse(UInt8, "00001000", base=2) # 0 -> point is not inside a bubble; 1 -> point is inside a bubble

abstract type AbstractPointTypes end

"""
    PointTypes{T <: AbstractFloat} <: AbstractPointTypes

Stores the point type information in a three dimensional cylindrical grid and also the position of the grid points.

# Fields
- `r::T`: The r-axis of the grid in `m`.
- `φ::T`: The φ-axis of the grid in `rad`.
- `z::T`: The r-axis of the grid in `m`.
- `pointtypes::Array{PointType, 3}`: The point type of each grid point.

"""
struct PointTypes{T <: AbstractFloat} <: AbstractPointTypes
    "Units: m"
    r::Array{T, 1}
    "Units: rad"
    φ::Array{T, 1}  # use radiance
    "Units: m"
    z::Array{T, 1}

    pointtypes::Array{PointType, 3}
end


function PointTypes(nt::NamedTuple)
    T = typeof(ustrip(nt.edges.r[1]))
    PointTypes(
        convert(Array{T}, ustrip.(uconvert.(u"m", nt.edges.r))),
        convert(Array{T}, ustrip.(uconvert.(u"rad", nt.edges.phi))),
        convert(Array{T}, ustrip.(uconvert.(u"m", nt.edges.z))),
        convert(Array{PointType}, ustrip.(uconvert.(NoUnits, nt.values)))
    )
end

Base.convert(T::Type{PointTypes}, x::NamedTuple) = T(x)


function NamedTuple(point_types::PointTypes)
    (
        values = point_types.pointtypes,
        edges = (
            r = point_types.r * u"m",
            phi = point_types.φ * u"rad",
            z = point_types.z * u"m"
        )
    )
end

Base.convert(T::Type{NamedTuple}, x::PointTypes) = T(x)



size(pts::PointTypes) = size(pts.pointtypes)
getindex(pts::PointTypes, i::Int) = getindex(pts.pointtypes, i)
getindex(pts::PointTypes, I::Vararg{Int, N}) where {N} = getindex(pts.pointtypes, I...)


function extent_2D_grid_to_3D(pts::PointTypes, n::Int)::PointTypes
    l = length(pts.φ)
    T::Type = eltype(pts.r)
    if l != 1 error("Grid is not of length 1 -> Not 2D.") end
    cyclic::T = 2π
    φ_step = cyclic / n
    φ::Array{T, 1} = [(i - 1) * φ_step for i in 1:n]
    pts_new = PointTypes{T}(pts.r, φ, pts.z, zeros(PointType, length(pts.r), length(φ), length(pts.z)) )
    for i in 1:n
        r = ((i - 1) * l + 1):(i * l)
        pts_new.pointtypes[:, r, :] = pts.pointtypes[:, :, :]
    end
    return pts_new
end
