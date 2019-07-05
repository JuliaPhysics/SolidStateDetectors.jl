
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
    return Int[searchsortednearest( grid[:r].ticks, point.r), searchsortednearest(grid[:φ].ticks, point.φ), searchsortednearest(grid[:z].ticks, point.z)]
end
function find_closest_gridpoint(point::CartesianPoint{T}, grid::CylindricalGrid{T}) where {T <: SSDFloat}
    find_closest_gridpoint(CylindricalPoint(point),grid)
end

function find_closest_gridpoint(point::CartesianPoint{T}, grid::CartesianGrid{T}) where {T <: SSDFloat}
    return Int[searchsortednearest( grid[:x].ticks, point.x), searchsortednearest(grid[:y].ticks, point.y), searchsortednearest(grid[:z].ticks, point.z)]
end
function find_closest_gridpoint(point::CylindricalPoint{T}, grid::CartesianGrid{T}) where {T <: SSDFloat}
    find_closest_gridpoint(CartesianPoint(point),grid)
end



