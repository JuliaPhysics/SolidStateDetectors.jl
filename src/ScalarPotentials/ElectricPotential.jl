"""
    struct ElectricPotential{T, N, S, AT} <: AbstractArray{T, N}
        
Electric potential of the simulation in units of volt (V).
        
## Parametric types 
* `T`: Element type of `data`.
* `N`: Dimension of the `grid` and `data` array.  
* `S`: Coordinate system (`Cartesian` or `Cylindrical`).
* `AT`: Axes type.  
        
## Fields
* `data::Array{T, N}`: Array containing the values of the electric potential at the discrete points of the `grid`.
* `grid::Grid{T, N, S, AT}`: [`Grid`](@ref) defining the discrete points for which the electric potential is determined.
"""
struct ElectricPotential{T, N, S, AT} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S, AT}
end

@inline size(epot::ElectricPotential{T, N, S}) where {T, N, S} = size(epot.data)
@inline length(epot::ElectricPotential{T, N, S}) where {T, N, S} = length(epot.data)
@inline getindex(epot::ElectricPotential{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(epot.data, I...)
@inline getindex(epot::ElectricPotential{T, N, S}, i::Int) where {T, N, S} = getindex(epot.data, i)
@inline getindex(epot::ElectricPotential{T, N, S}, s::Symbol) where {T, N, S} = getindex(epot.grid, s)


function Base.NamedTuple(epot::ElectricPotential{T, 3}) where {T}
    return (
        grid = NamedTuple(epot.grid),
        values = epot.data * u"V",
    )
end
Base.convert(T::Type{NamedTuple}, x::ElectricPotential) = T(x)
    
function ElectricPotential(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(ustrip(nt.values[1]))
    S = get_coordinate_system(grid)
    N = get_number_of_dimensions(grid)
    ElectricPotential{T, N, S, typeof(grid.axes)}( ustrip.(uconvert.(u"V", nt.values)), grid)
end
Base.convert(T::Type{ElectricPotential}, x::NamedTuple) = T(x)

"""
    get_ticks_at_positions_of_large_gradient(epot::ElectricPotential)

The electric potential is analyzed in order to find and return ticks where the gradient (electric field) is strong (relativ to its maximum).

In the 1D case of a pn-junction, the electric field strength is the largest at the position of the pn-junction.
Thus, this function is likely to return ticks which are located close to the pn-junction of a semiconductor.
"""
function get_ticks_at_positions_of_large_gradient(epot::ElectricPotential{T, 3}) where T
    dims = (1, 2, 3)
    mps::NTuple{3, Vector{T}} = broadcast(idim -> StatsBase.midpoints(epot.grid[idim]), dims)
    extended = length.(mps) .> 0
    max_diffs = broadcast(idim -> extended[idim] ? map(i->maximum(abs.(selectdim(epot.data, idim, i+1) .- selectdim(epot.data, idim, i))), eachindex(mps[idim]))::Vector{T} : T[], dims);
    max_diff = broadcast(idim -> extended[idim] ? maximum(max_diffs[idim]) : T(0), dims)
    normalized_max_diffs = broadcast(idim -> extended[idim] ? max_diffs[idim] ./ max_diff[idim] : T[], dims)
    inds = broadcast(idim -> extended[idim] ? findall(Δ -> Δ > 0.5, normalized_max_diffs[idim]) : Int[], dims) # 0.5 seems to be a good limit as we don't want to add too many ticks. 
    return broadcast(idim -> extended[idim] ? mps[idim][inds[idim]] : T[], dims)
end