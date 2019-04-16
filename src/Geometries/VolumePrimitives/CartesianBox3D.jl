"""
    mutable struct CartesianBox3D{T} <: AbstractGeometry{T, 3}

Very simple rectengular box in cartesian coordinates.
"""
struct CartesianBox3D{T} <: AbstractGeometry{T, 3} # should be immutable, mutable right now just for testing.
    hierarchy::T
    x::Tuple{T, T}
    y::Tuple{T, T}
    z::Tuple{T, T}
end

function in(pt::CartesianPoint{T}, g::CartesianBox3D{T})::Bool where {T}
    return (g.x[1] <= pt[1] <= g.x[2]) && (g.y[1] <= pt[2] <= g.y[2]) && (g.z[1] <= pt[3] <= g.z[2])
end

@inline in(pt::CylindricalPoint, g::CartesianBox3D)::Bool = in(CartesianPoint(pt), g)

function CartesianBox3D{T}(dict::Dict{Any, Any}, inputunit::Unitful.Units)::CartesianBox3D{T} where {T <: SSDFloat}
    haskey(dict, "hierarchy") ? h::Int = dict["hierarchy"] : h = 1
    return CartesianBox3D{T}(
        h,
        (geom_round(ustrip(uconvert(u"m", T(dict["xStart"]) * inputunit ))), geom_round(ustrip(uconvert(u"m", T(dict["xStop"]) * inputunit)))),
        (geom_round(ustrip(uconvert(u"m", T(dict["yStart"]) * inputunit ))), geom_round(ustrip(uconvert(u"m", T(dict["yStop"]) * inputunit)))),
        (geom_round(ustrip(uconvert(u"m", T(dict["zStart"]) * inputunit ))), geom_round(ustrip(uconvert(u"m", T(dict["zStop"]) * inputunit))))  )
end

function Geometry(T::DataType, t::Val{:CartesianBox3D}, dict::Dict{Any, Any}, inputunit::Unitful.Units)
    return CartesianBox3D{T}(dict, inputunit)
end

function get_important_points(g::CartesianBox3D{T})::NTuple{3, Vector{T}} where {T <: SSDFloat}
    v1::Vector{T} = T[g.x[1], g.x[2]] #[g.x[1], g.x[2]]
    v2::Vector{T} = T[g.y[1], g.y[2]]
    v3::Vector{T} = T[g.z[1], g.z[2]]
    return v1, v2, v3
end

function get_important_points(g::CartesianBox3D{T}, ::Val{:r})::Vector{T} where {T <: SSDFloat}
    @warn "Not yet implemented"
    return T[]
end
function get_important_points(g::CartesianBox3D{T}, ::Val{:φ})::Vector{T} where {T <: SSDFloat}
    @warn "Not yet implemented"
    return T[]
end
function get_important_points(g::CartesianBox3D{T}, ::Val{:z})::Vector{T} where {T <: SSDFloat}
    return T[g.z[1], g.z[2]]
end
function get_important_points(g::CartesianBox3D{T}, ::Val{:x})::Vector{T} where {T <: SSDFloat}
    return T[g.x[1], g.x[2]]
end
function get_important_points(g::CartesianBox3D{T}, ::Val{:y})::Vector{T} where {T <: SSDFloat}
    return T[g.y[1], g.y[2]]
end


function vertices(cb::CartesianBox3D{T})::Vector{CartesianPoint{T}} where {T <: SSDFloat}
    v::Vector{CartesianPoint{T}} = CartesianPoint{T}[
        CartesianPoint{T}(cb.x[1], cb.y[1], cb.z[1]),
        CartesianPoint{T}(cb.x[2], cb.y[1], cb.z[1]),
        CartesianPoint{T}(cb.x[2], cb.y[2], cb.z[1]),
        CartesianPoint{T}(cb.x[1], cb.y[2], cb.z[1]),
        CartesianPoint{T}(cb.x[1], cb.y[1], cb.z[2]),
        CartesianPoint{T}(cb.x[2], cb.y[1], cb.z[2]),
        CartesianPoint{T}(cb.x[2], cb.y[2], cb.z[2]),
        CartesianPoint{T}(cb.x[1], cb.y[2], cb.z[2]),
    ]
    return v
end

function LineSegments(cb::CartesianBox3D{T})::Vector{LineSegment{T, 3, :Cartesian}} where {T}
    v::Vector{CartesianPoint{T}} = vertices(cb)
    return LineSegment{T, 3, :Cartesian}[
        LineSegment(v[1], v[2]),
        LineSegment(v[2], v[3]),
        LineSegment(v[3], v[4]),
        LineSegment(v[4], v[1]),
        LineSegment(v[1], v[5]),
        LineSegment(v[2], v[6]),
        LineSegment(v[3], v[7]),
        LineSegment(v[4], v[8]),
        LineSegment(v[5], v[6]),
        LineSegment(v[6], v[7]),
        LineSegment(v[7], v[8]),
        LineSegment(v[8], v[5])
    ]
end

@recipe function f(cb::CartesianBox3D{T}) where {T <: SSDFloat}
    ls = LineSegments(cb)
    @series begin
        ls
    end
end

function sample(cb::CartesianBox3D{T}, stepsize::Vector{T})  where {T <: SSDFloat}
    samples  = CartesianPoint{T}[]
    for x in cb.x[1] : stepsize[1] : cb.x[2]
        for y in cb.y[1] :stepsize[2] : cb.y[2]
            for z in cb.z[1] : stepsize[3] : cb.z[2]
                push!(samples,CartesianPoint{T}(x, y, z))
            end
        end
    end
    return samples
end
