abstract type AbstractGrid{T, N} <: AbstractArray{T, N} end

"""
    T: tick type
    N: N dimensional
    S: System (Cartesian, Cylindrical...)
"""
struct Grid{T, N, S} <: AbstractGrid{T, N}
    axes::NTuple{N, DiscreteAxis{T, BL, BR, I} where {BL, BR, I}} 
end

const CartesianGrid{T, N} = Grid{T, N, :cartesian} 
const CartesianGrid1D{T} = CartesianGrid{T, 1}
const CartesianGrid2D{T} = CartesianGrid{T, 2}
const CartesianGrid3D{T} = CartesianGrid{T, 3}
const RadialGrid{T} = Grid{T, 1, :Radial} 
const PolarGrid{T} = Grid{T, 2, :Polar} 
const CylindricalGrid{T} = Grid{T, 3, :cylindrical} 
const SphericalGrid{T} = Grid{T, 3, :Spherical} 

@inline size(g::Grid{T, N, S}) where {T, N, S} = size.(g.axes, 1)
@inline length(g::Grid{T, N, S}) where {T, N, S} = prod(size(g))
@inline getindex(g::Grid{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = broadcast(getindex, g.axes, I)
@inline getindex(g::Grid{T, N, S}, i::Int) where {T, N, S} = getproperty(g, :axes)[i]
@inline getindex(g::Grid{T, N, S}, s::Symbol) where {T, N, S} = getindex(g, Val{s}())

@inline getproperty(g::Grid{T, N, S}, s::Symbol) where {T, N, S} = getproperty(g, Val{s}())
@inline getproperty(g::Grid{T}, ::Val{:axes}) where {T} = getfield(g, :axes)

@inline getproperty(g::CylindricalGrid{T}, ::Val{:axes}) where {T} = getfield(g, :axes)
@inline getproperty(g::CylindricalGrid{T}, ::Val{:r}) where {T} = @inbounds g.axes[1]
@inline getproperty(g::CylindricalGrid{T}, ::Val{:φ}) where {T} = @inbounds g.axes[2]
@inline getproperty(g::CylindricalGrid{T}, ::Val{:z}) where {T} = @inbounds g.axes[3]
@inline getproperty(g::CartesianGrid{T}, ::Val{:x}) where {T} = @inbounds g.axes[1]
@inline getproperty(g::CartesianGrid{T}, ::Val{:y}) where {T} = @inbounds g.axes[2]
@inline getproperty(g::CartesianGrid{T}, ::Val{:z}) where {T} = @inbounds g.axes[3]

@inline getindex(g::CylindricalGrid{T}, ::Val{:r}) where {T} = @inbounds g.axes[1]
@inline getindex(g::CylindricalGrid{T}, ::Val{:φ}) where {T} = @inbounds g.axes[2]
@inline getindex(g::CylindricalGrid{T}, ::Val{:z}) where {T} = @inbounds g.axes[3]
@inline getindex(g::CartesianGrid{T}, ::Val{:x}) where {T} = @inbounds g.axes[1]
@inline getindex(g::CartesianGrid{T}, ::Val{:y}) where {T} = @inbounds g.axes[2]
@inline getindex(g::CartesianGrid{T}, ::Val{:z}) where {T} = @inbounds g.axes[3]

function sizeof(g::Grid{T, N, S}) where {T, N, S}
    return sum( sizeof.(g.axes) )
end

function print(io::IO, g::Grid{T, N, S}) where {T, N, S}
    print(io, "Grid{$T, $N, $S}", g.axes)
end
function println(io::IO, g::Grid{T, N, S}) where {T, N, S}
    println(" Grid{$T, $N, $S}")
    for (i, ax) in enumerate(g.axes)
        println(io, "  Axis $(i): ", ax)
    end
end
show(io::IO, g::Grid{T, N, S}) where {T, N, S} = print(io, g)
show(io::IO, ::MIME"text/plain", g::Grid{T, N, S}) where {T, N, S} = show(io, g)


function searchsortednearest(pt::CylindricalPoint{T}, g::CylindricalGrid{T})::NTuple{3, Int} where {T <: SSDFloat}
    ir::Int = searchsortednearest(g.axes[1].ticks, pt.r)
    iφ::Int = searchsortednearest(g.axes[2].ticks, pt.φ)
    iz::Int = searchsortednearest(g.axes[3].ticks, pt.z)
    ir, iφ, iz
end
function searchsortednearest(pt::CartesianPoint{T}, g::CylindricalGrid{T})::NTuple{3, Int} where {T <: SSDFloat}
    return searchsortednearest(CylindricalPoint(pt), g)
end

function searchsortednearest(pt::CartesianPoint{T}, g::CartesianGrid{T, 3})::NTuple{3, Int} where {T <: SSDFloat}
    ix::Int = searchsortednearest(g.axes[1].ticks, pt.x)
    iy::Int = searchsortednearest(g.axes[2].ticks, pt.y)
    iz::Int = searchsortednearest(g.axes[3].ticks, pt.z)
    ix, iy, iz
end
function searchsortednearest(pt::CylindricalPoint{T}, g::CartesianGrid{T, 3})::NTuple{3, Int} where {T <: SSDFloat}
    return searchsortednearest(CartesianPoint(pt), g)
end


function check_grid(grid::CylindricalGrid{T})::Nothing where {T}
    nr::Int, nφ::Int, nz::Int = size(grid)
    @assert iseven(nz) "GridError: Field simulation algorithm in cylindrical coordinates need an even number of grid points in z. This is not the case. #z-ticks = $(nz)."
    @assert (iseven(nφ) || (nφ == 1)) "GridError: Field simulation algorithm in cylindrical coordinates need an even number of grid points in φ or just one point (2D). This is not the case. #φ-ticks = $(nφ)."
    return nothing
end

function check_grid(grid::CartesianGrid3D{T})::Nothing where {T}
    nx::Int, ny::Int, nz::Int = size(grid)
    @assert iseven(nx) "GridError: Field simulation algorithm in cartesian coordinates need an even number of grid points in x. This is not the case. #x-ticks = $(nx)."
    return nothing
end

include("RefineGrid.jl")

@recipe function f(grid::Grid{T, N, S}) where {T, N, S}
    layout --> N
    st --> :stephist
    legend --> false

    for iax in 1:N
        @series begin
            subplot := iax
            nbins --> div(length(grid[iax]), 2)
            xlabel --> "Grid point density - Axis $(iax)"
            grid[iax]
        end
    end
end


function get_coordinate_system(grid::Grid{T, N, S}) where {T, N, S}
    return S
end
function get_number_of_dimensions(grid::Grid{T, N, S}) where {T, N, S}
    return N
end
function eltype(grid::Grid{T, N, S})::DataType where {T, N, S}
    return T
end

function get_boundary_types(grid::Grid{T, N, S}) where {T, N, S}
   return get_boundary_types.(grid.axes) 
end


function Grid(nt::NamedTuple)
    if nt.coordtype == "cylindrical"
        axr::DiscreteAxis = DiscreteAxis(nt.axes.r, unit=u"m")
        axφ::DiscreteAxis = DiscreteAxis(nt.axes.phi, unit=u"rad")
        axz::DiscreteAxis = DiscreteAxis(nt.axes.z, unit=u"m")
        T = typeof(axr.ticks[1])
        return Grid{T, 3, :cylindrical}( (axr, axφ, axz) )
    elseif nt.coordtype == "cartesian"
        axx::DiscreteAxis = DiscreteAxis(nt.axes.x, unit=u"m")
        axy::DiscreteAxis = DiscreteAxis(nt.axes.y, unit=u"m")
        axz = DiscreteAxis(nt.axes.z, unit=u"m")
        T = typeof(axx.ticks[1])
        return Grid{T, 3, :cartesian}( (axx, axy, axz) )
    else
        error("`coordtype` = $(nt.coordtype) is not valid.")
    end
end

Base.convert(T::Type{Grid}, x::NamedTuple) = T(x)

function NamedTuple(grid::Grid{T, 3, :cylindrical}) where {T}
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
function NamedTuple(grid::Grid{T, 3, :cartesian}) where {T}
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

