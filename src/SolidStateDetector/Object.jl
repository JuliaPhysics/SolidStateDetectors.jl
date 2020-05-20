
abstract type AbstractObject{T <: SSDFloat} end


@inline function in(pt::AbstractCoordinatePoint{T}, c::AbstractObject{T})::Bool where {T}
    return in(pt, c.geometry)
end

# function in(p::AbstractCoordinatePoint{T}, c::AbstractObject{T}, rs::Vector{T}) where {T <: SSDFloat}
#     rv = false
#     for g in c.geometry_positive
#         if typeof(g) in []#[ConeMantle{T}]
#             in(p, g, rs) ? rv = true : nothing
#         else
#             (p in g) ? rv = true : nothing
#         end
#     end
#     return rv
# end

function in(p::AbstractCoordinatePoint{T}, v::AbstractVector{<:AbstractObject{T}}) where {T <: SSDFloat}
    rv = false
    for object in v
        if p in object
            rv = true
        end
    end
    return rv
end

function find_closest_gridpoint(point::CylindricalPoint{T}, grid::CylindricalGrid{T}) where {T <: SSDFloat}
    return Int[searchsortednearest( grid.axes[1].ticks, point.r), searchsortednearest(grid.axes[2].ticks, point.φ), searchsortednearest(grid.axes[3].ticks, point.z)]
end
function find_closest_gridpoint(point::CartesianPoint{T}, grid::CylindricalGrid{T}) where {T <: SSDFloat}
    find_closest_gridpoint(CylindricalPoint(point),grid)
end

function find_closest_gridpoint(point::CartesianPoint{T}, grid::CartesianGrid{T}) where {T <: SSDFloat}
    @inbounds return Int[searchsortednearest( grid.axes[1].ticks, point.x), searchsortednearest(grid.axes[2].ticks, point.y), searchsortednearest(grid.axes[3].ticks, point.z)]
end
function find_closest_gridpoint(point::CylindricalPoint{T}, grid::CartesianGrid{T}) where {T <: SSDFloat}
    find_closest_gridpoint(CartesianPoint(point),grid)
end



