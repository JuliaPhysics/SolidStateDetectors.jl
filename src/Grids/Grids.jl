abstract type AbstractGrid{T, N} <: AbstractArray{T, N} end

"""
    struct Grid{T, N, S <: AbstractCoordinateSystem, AT} <: AbstractGrid{T, N}

Collection of `N` axes that define the dimensions of the grid needed to calculate 
[`ElectricPotential`](@ref), [`ElectricField`](@ref) or [`WeightingPotential`](@ref).

## Parametric types
* `T`: Tick type (element type) of the axes.
* `N`: Dimension of the grid.
* `S`: Coordinate system (`Cartesian` or `Cylindrical`).
* `AT`: Axes type.

## Fields
* `axes::AT`: Tuple of length `N` containing `DiscreteAxis` for each dimension of the grid.

See also [`DiscreteAxis`](@ref).
"""
struct Grid{T, N, S <: AbstractCoordinateSystem, AT} <: AbstractGrid{T, N}
    axes::AT
end

const CartesianGrid{T, N} = Grid{T, N, Cartesian}
const CartesianGrid1D{T} = CartesianGrid{T, 1}
const CartesianGrid2D{T} = CartesianGrid{T, 2}
const CartesianGrid3D{T} = CartesianGrid{T, 3}
const CylindricalGrid{T} = Grid{T, 3, Cylindrical}
#const RadialGrid{T} = Grid{T, 1, Radial}
#const PolarGrid{T} = Grid{T, 2, Polar}
#const SphericalGrid{T} = Grid{T, 3, Spherical}

CylindricalGrid{T}(a) where {T} = Grid{T, 3, Cylindrical, typeof(a)}(a)
CartesianGrid3D{T}(a) where {T} = Grid{T, 3, Cartesian, typeof(a)}(a)

@inline size(grid::Grid{T, N, S}) where {T, N, S} = size.(grid.axes, 1)
@inline length(grid::Grid{T, N, S}) where {T, N, S} = prod(size(grid))
@inline getindex(grid::Grid{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = broadcast(getindex, grid.axes, I)
@inline getindex(grid::Grid{T, N, S}, i::Int) where {T, N, S} = getproperty(grid, :axes)[i]
@inline getindex(grid::Grid{T, N, S}, s::Symbol) where {T, N, S} = getindex(grid, Val{s}())

@inline getproperty(grid::Grid{T, N, S}, s::Symbol) where {T, N, S} = getproperty(grid, Val{s}())
@inline getproperty(grid::Grid{T}, ::Val{:axes}) where {T} = getfield(grid, :axes)

@inline getproperty(grid::CylindricalGrid{T}, ::Val{:axes}) where {T} = getfield(grid, :axes)
@inline getproperty(grid::CylindricalGrid{T}, ::Val{:r}) where {T} = @inbounds grid.axes[1]
@inline getproperty(grid::CylindricalGrid{T}, ::Val{:φ}) where {T} = @inbounds grid.axes[2]
@inline getproperty(grid::CylindricalGrid{T}, ::Val{:z}) where {T} = @inbounds grid.axes[3]
@inline getproperty(grid::CartesianGrid{T}, ::Val{:x}) where {T} = @inbounds grid.axes[1]
@inline getproperty(grid::CartesianGrid{T}, ::Val{:y}) where {T} = @inbounds grid.axes[2]
@inline getproperty(grid::CartesianGrid{T}, ::Val{:z}) where {T} = @inbounds grid.axes[3]

@inline getindex(grid::CylindricalGrid{T}, ::Val{:r}) where {T} = @inbounds grid.axes[1]
@inline getindex(grid::CylindricalGrid{T}, ::Val{:φ}) where {T} = @inbounds grid.axes[2]
@inline getindex(grid::CylindricalGrid{T}, ::Val{:z}) where {T} = @inbounds grid.axes[3]
@inline getindex(grid::CartesianGrid{T}, ::Val{:x}) where {T} = @inbounds grid.axes[1]
@inline getindex(grid::CartesianGrid{T}, ::Val{:y}) where {T} = @inbounds grid.axes[2]
@inline getindex(grid::CartesianGrid{T}, ::Val{:z}) where {T} = @inbounds grid.axes[3]

function sizeof(grid::Grid{T, N, S}) where {T, N, S}
    return sum( sizeof.(grid.axes) )
end

function print(io::IO, grid::Grid{T, N, S}) where {T, N, S}
    print(io, "Grid{$T, $N, $S}", grid.axes)
end
function println(io::IO, grid::Grid{T, N, S}) where {T, N, S}
    println(" Grid{$T, $N, $S}")
    for (i, ax) in enumerate(grid.axes)
        println(io, "  Axis $(i): ", ax)
    end
end
show(io::IO, grid::Grid{T, N, S}) where {T, N, S} = print(io, grid)
show(io::IO, ::MIME"text/plain", grid::Grid{T, N, S}) where {T, N, S} = show(io, grid)

function check_grid(grid::CylindricalGrid{T})::Nothing where {T}
    nr::Int, nφ::Int, nz::Int = size(grid)
    @assert iseven(nz) "GridError: Field simulation algorithm in cylindrical coordinates needs an even number of grid points in z. This is not the case. #z-ticks = $(nz)."
    @assert (iseven(nφ) || (nφ == 1)) "GridError: Field simulation algorithm in cylindrical coordinates needs an even number of grid points in φ or just one point (2D). This is not the case. #φ-ticks = $(nφ)."
    return nothing
end

function check_grid(grid::CartesianGrid3D{T})::Nothing where {T}
    nx::Int, ny::Int, nz::Int = size(grid)
    @assert iseven(nx) "GridError: Field simulation algorithm in cartesian coordinates needs an even number of grid points in x. This is not the case. #x-ticks = $(nx)."
    return nothing
end

function get_coordinate_system(grid::Grid{T, N, S}) where {T, N, S}
    return S
end
function get_number_of_dimensions(grid::Grid{T, N, S}) where {T, N, S}
    return N
end
function Base.eltype(grid::Grid{T, N, S})::DataType where {T, N, S}
    return T
end

function get_boundary_types(grid::Grid{T, N, S}) where {T, N, S}
   return get_boundary_types.(grid.axes)
end

TicksTuple(grid::Grid{T, 3, Cartesian}) where {T} = (x = grid.axes[1].ticks, y = grid.axes[2].ticks, z = grid.axes[3].ticks)
TicksTuple(grid::Grid{T, 3, Cylindrical}) where {T} = (r = grid.axes[1].ticks, φ = grid.axes[2].ticks, z = grid.axes[3].ticks)

function Grid(nt::NamedTuple)
    if nt.coordtype == "cylindrical"
        axr::DiscreteAxis = DiscreteAxis(nt.axes.r, unit=u"m")
        axφ::DiscreteAxis = DiscreteAxis(nt.axes.phi, unit=u"rad")
        axz::DiscreteAxis = DiscreteAxis(nt.axes.z, unit=u"m")
        T = typeof(axr.ticks[1])
        return Grid{T, 3, Cylindrical}( (axr, axφ, axz) )
    elseif nt.coordtype == "cartesian"
        axx::DiscreteAxis = DiscreteAxis(nt.axes.x, unit=u"m")
        axy::DiscreteAxis = DiscreteAxis(nt.axes.y, unit=u"m")
        axz = DiscreteAxis(nt.axes.z, unit=u"m")
        T = typeof(axx.ticks[1])
        return Grid{T, 3, Cartesian}( (axx, axy, axz) )
    else
        error("`coordtype` = $(nt.coordtype) is not valid.")
    end
end

Base.convert(T::Type{Grid}, x::NamedTuple) = T(x)

function NamedTuple(grid::Grid{T, 3, Cylindrical}) where {T}
    axr::DiscreteAxis{T} = grid.axes[1]
    axφ::DiscreteAxis{T} = grid.axes[2]
    axz::DiscreteAxis{T} = grid.axes[3]
    return (
        coordtype = "cylindrical",
        ndims = 3,
        axes = (
            r   = NamedTuple(axr, unit=u"m"),
            phi = NamedTuple(axφ, unit=u"rad"),
            z   = NamedTuple(axz, unit=u"m"),
        )
    )
end
function NamedTuple(grid::Grid{T, 3, Cartesian}) where {T}
    axx::DiscreteAxis{T} = grid.axes[1]
    axy::DiscreteAxis{T} = grid.axes[2]
    axz::DiscreteAxis{T} = grid.axes[3]
    return (
        coordtype = "cartesian",
        ndims = 3,
        axes = (
            x = NamedTuple(axx, unit=u"m"),
            y = NamedTuple(axy, unit=u"m"),
            z = NamedTuple(axz, unit=u"m"),
        )
    )
end

Base.convert(T::Type{NamedTuple}, x::Grid) = T(x)


function find_closest_gridpoint(pt::CylindricalPoint{T}, grid::CylindricalGrid{T})::NTuple{3,Int} where {T <: SSDFloat}
    return (searchsortednearest(grid.axes[1].ticks, pt.r), searchsortednearest(grid.axes[2].ticks, pt.φ), searchsortednearest(grid.axes[3].ticks, pt.z))
end
function find_closest_gridpoint(pt::CartesianPoint{T}, grid::CylindricalGrid{T})::NTuple{3,Int} where {T <: SSDFloat}
    find_closest_gridpt(CylindricalPoint(pt),grid)
end

function find_closest_gridpoint(pt::CartesianPoint{T}, grid::CartesianGrid{T})::NTuple{3,Int} where {T <: SSDFloat}
    @inbounds return (searchsortednearest(grid.axes[1].ticks, pt.x), searchsortednearest(grid.axes[2].ticks, pt.y), searchsortednearest(grid.axes[3].ticks, pt.z))
end
function find_closest_gridpoint(pt::CylindricalPoint{T}, grid::CartesianGrid{T})::NTuple{3,Int} where {T <: SSDFloat}
    find_closest_gridpt(CartesianPoint(pt),grid)
end
