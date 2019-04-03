"""
    mutable struct CartesianBox3D{T} <: AbstractGeometry{T, 3}

Very simple rectengular box in cartesian coordinates.
"""
mutable struct CartesianBox3D{T} <: AbstractGeometry{T, 3} # should be immutable, mutable right now just for testing.
    x::Tuple{T, T}
    y::Tuple{T, T}
    z::Tuple{T, T}
end

function in(pt::CartesianPoint{T}, g::CartesianBox3D{T})::Bool where {T}
    return (g.x[1] <= pt[1] <= g.x[2]) && (g.y[1] <= pt[2] <= g.y[2]) && (g.z[1] <= pt[3] <= g.z[2])
end

@inline in(pt::CylindricalPoint, g::CartesianBox3D)::Bool = in(CartesianPoint(pt), g)

function CartesianBox3D{T}(dict::Dict{Any, Any}, inputunit::Unitful.Units)::CartesianBox3D{T} where {T <: AbstractFloat}
    return CartesianBox3D{T}(
        (geom_round(ustrip(uconvert(u"m", T(dict["xStart"]) * inputunit ))), geom_round(ustrip(uconvert(u"m", T(dict["xStop"]) * inputunit)))),
        (geom_round(ustrip(uconvert(u"m", T(dict["yStart"]) * inputunit ))), geom_round(ustrip(uconvert(u"m", T(dict["yStop"]) * inputunit)))),
        (geom_round(ustrip(uconvert(u"m", T(dict["zStart"]) * inputunit ))), geom_round(ustrip(uconvert(u"m", T(dict["zStop"]) * inputunit))))  )
end

function Geometry(T::DataType, t::Val{:CartesianBox3D}, dict::Dict{Any, Any}, inputunit::Unitful.Units)
    return CartesianBox3D{T}(dict, inputunit)
end

function get_important_points(g::CartesianBox3D{T})::NTuple{3, Vector{T}} where {T <: AbstractFloat}
    v1::Vector{T} = T[g.x[1], g.x[2]] #[g.x[1], g.x[2]]
    v2::Vector{T} = T[g.y[1], g.y[2]]
    v3::Vector{T} = T[g.z[1], g.z[2]]
    return v1, v2, v3
end

function get_important_points(g::CartesianBox3D{T}, ::Val{:r})::Vector{T} where {T <: AbstractFloat}
    @warn "Not yet implemented"
    return T[]
end
function get_important_points(g::CartesianBox3D{T}, ::Val{:φ})::Vector{T} where {T <: AbstractFloat}
    @warn "Not yet implemented"
    return T[]
end
function get_important_points(g::CartesianBox3D{T}, ::Val{:z})::Vector{T} where {T <: AbstractFloat}
    return T[g.z[1], g.z[2]]
end
function get_important_points(g::CartesianBox3D{T}, ::Val{:x})::Vector{T} where {T <: AbstractFloat}
    return T[g.x[1], g.x[2]]
end
function get_important_points(g::CartesianBox3D{T}, ::Val{:y})::Vector{T} where {T <: AbstractFloat}
    return T[g.y[1], g.y[2]]
end

