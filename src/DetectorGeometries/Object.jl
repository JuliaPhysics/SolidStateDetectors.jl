
abstract type AbstractObject{T <: SSDFloat} end


@inline function in(pt::AbstractCoordinatePoint{T}, c::AbstractObject{T})::Bool where {T}
    return in(pt, c.geometry)
end

# function in(p::AbstractCoordinatePoint{T}, c::AbstractObject{T}, rs::Vector{T}) where T
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

function in(p::AbstractCoordinatePoint{T}, v::AbstractVector{<:AbstractObject{T}}) where T
    rv = false
    for object in v
        if p in object
            rv = true
        end
    end
    return rv
end

function find_closest_gridpoint(point::CylindricalPoint{T}, grid::CylindricalGrid{T}) where T
    return Int[searchsortednearest( grid[:r].ticks, point.r), searchsortednearest(grid[:φ].ticks, point.φ), searchsortednearest(grid[:z].ticks, point.z)]
end
function find_closest_gridpoint(point::CartesianPoint{T}, grid::CylindricalGrid{T}) where T
    find_closest_gridpoint(CylindricalPoint(point),grid)
end

function find_closest_gridpoint(point::CartesianPoint{T}, grid::CartesianGrid{T}) where T
    return Int[searchsortednearest( grid[:x].ticks, point.x), searchsortednearest(grid[:y].ticks, point.y), searchsortednearest(grid[:z].ticks, point.z)]
end
function find_closest_gridpoint(point::CylindricalPoint{T}, grid::CartesianGrid{T}) where T
    find_closest_gridpoint(CartesianPoint(point),grid)
end

function paint_object(object::AbstractObject, grid::CylindricalGrid{T}) where T
    stepsize::Vector{T}= [minimum(diff(grid[:r].ticks)), IntervalSets.width(grid[:φ].interval) == 0.0 ? 0.05236 : minimum(diff(grid[:φ].ticks)), minimum(diff(grid[:z].ticks))]
    samples = filter(x-> x in object.geometry, vcat([sample(g, stepsize) for g in object.geometry_positive]...))
    object_gridpoints = unique!([find_closest_gridpoint(sample_point,grid) for sample_point in samples])
    return object_gridpoints
end
function paint_object(object::AbstractObject{T}, grid::CylindricalGrid{T}, φ::T) where T
    closest_φ_idx=searchsortednearest(grid[:φ].ticks, φ)
    stepsize::Vector{T}= [minimum(diff(grid[:r].ticks)), IntervalSets.width(grid[:φ].interval) == 0.0 ? 0.05236 : minimum(diff(grid[:φ].ticks)), minimum(diff(grid[:z].ticks))]
    samples = filter(x-> x in object.geometry, vcat([sample(g, stepsize) for g in object.geometry_positive]...))
    object_gridpoints = unique!([find_closest_gridpoint(sample_point,grid) for sample_point in samples])
    return filter(x -> x[2]==closest_φ_idx, object_gridpoints)
end

function paint_object(object::AbstractObject{T}, grid::CartesianGrid{T}) where T
    stepsize::Vector{T}= [minimum(diff(grid[:x].ticks)), minimum(diff(grid[:y].ticks)), minimum(diff(grid[:z].ticks))]
    samples = filter(x-> x in object.geometry, vcat([sample(g, stepsize) for g in object.geometry_positive]...))
    object_gridpoints = unique!([find_closest_gridpoint(sample_point,grid) for sample_point in samples])
    return object_gridpoints
end
