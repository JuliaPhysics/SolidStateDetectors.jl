"""
    struct ImpurityScale{T, N, S, AT} <: AbstractArray{T, N}
        
Impurity scalar field of the simulation.
It is kinda an alpha map for the impurity density of semiconductors.
It can takes values between `0` and `1`:
 * `1`: The impurity density has its full value. For grid points in depleted region of the semiconductor. 
 * `]0,1[`: The impurity density is scaled down but not zero. For grid points at the edge of the depleted region which are partially depleted.
 * `0`: The impurity density is set to `0`. For grid points in undepleted regions of the semiconductor.

        
## Parametric types 
* `T`: Element type of `data`.
* `N`: Dimension of the `grid` and `data` array.  
* `S`: Coordinate system (`Cartesian` or `Cylindrical`).
* `AT`: Axes type.  
        
## Fields
* `data::Array{T, N}`: Array containing the values of the impurity scale at the discrete points of the `grid`.
* `grid::Grid{T, N, S, AT}`: [`Grid`](@ref) defining the discrete points for which the impurity scale is determined.
"""
struct ImpurityScale{T, N, S, AT} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S, AT}
end

@inline size(impscale::ImpurityScale{T, N, S}) where {T, N, S} = size(impscale.data)
@inline length(impscale::ImpurityScale{T, N, S}) where {T, N, S} = length(impscale.data)
@inline getindex(impscale::ImpurityScale{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(impscale.data, I...)
@inline getindex(impscale::ImpurityScale{T, N, S}, i::Int) where {T, N, S} = getindex(impscale.data, i)
@inline getindex(impscale::ImpurityScale{T, N, S}, s::Symbol) where {T, N, S} = getindex(impscale.grid, s)


function Base.NamedTuple(impscale::ImpurityScale{T, 3}) where {T}
    return (
        grid = NamedTuple(impscale.grid),
        values = impscale.data,
    )
end
Base.convert(T::Type{NamedTuple}, x::ImpurityScale) = T(x)
    
function ImpurityScale(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(ustrip(nt.values[1]))
    S = get_coordinate_system(grid)
    N = get_number_of_dimensions(grid)
    ImpurityScale{T, N, S, typeof(grid.axes)}(nt.values, grid)
end
Base.convert(T::Type{ImpurityScale}, x::NamedTuple) = T(x)


"""
    get_ticks_at_positions_of_edge_of_depleted_volumes(impscale::ImpurityScale{T, 3})

The impurity scale field is analyzed in order to find and return ticks where the gradient is strong,
which is the case at the surface of the depleted volume of a semiconductor.
"""
function get_ticks_at_positions_of_edge_of_depleted_volumes(impscale::ImpurityScale{T, 3}) where T
    dims = (1, 2, 3)
    mps::NTuple{3, Vector{T}} = broadcast(idim -> StatsBase.midpoints(impscale.grid[idim]), dims)
    extended = length.(mps) .> 0
    max_diffs = broadcast(idim -> extended[idim] ? map(i->maximum(abs.(selectdim(impscale.data, idim, i+1) .- selectdim(impscale.data, idim, i))), eachindex(mps[idim]))::Vector{T} : T[], dims);
    inds = broadcast(idim -> extended[idim] ? findall(Δ -> 0.5 < Δ < 1, max_diffs[idim]) : Int[], dims) # 0.5 seems to be a good limit as we don't want to add too many ticks. 
    # As the impurity scale only has values between `0` and `1` the change (`max_diffs`) can also only between `0` and `1`.
    # We want to exclude a change of `1` as this happens at the surface of a semiconductor where the detector is undepleted
    # and we only want to generate ticks in the semiconductor at the position between undepleted and depleted regions.
    return broadcast(idim -> extended[idim] ? mps[idim][inds[idim]] : T[], dims)
end
