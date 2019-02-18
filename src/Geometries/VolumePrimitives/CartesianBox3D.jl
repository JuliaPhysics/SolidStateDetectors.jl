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

function CartesianBox3D{T}(dict::Dict{String, Any}, inputunit::Unitful.Units)::CartesianBox3D{T} where {T <: AbstractFloat}
    return CartesianBox3D{T}(
        (geom_round(ustrip(uconvert(u"m", T(dict["xStart"]) * inputunit ))), geom_round(ustrip(uconvert(u"m", T(dict["xStop"]) * inputunit)))),
        (geom_round(ustrip(uconvert(u"m", T(dict["yStart"]) * inputunit ))), geom_round(ustrip(uconvert(u"m", T(dict["yStop"]) * inputunit)))),
        (geom_round(ustrip(uconvert(u"m", T(dict["zStart"]) * inputunit ))), geom_round(ustrip(uconvert(u"m", T(dict["zStop"]) * inputunit))))  )
end

function Geometry(T::DataType, t::Val{:CartesianBox3D}, dict::Dict{String, Any}, inputunit::Unitful.Units)
    return CartesianBox3D{T}(dict, inputunit)
end

