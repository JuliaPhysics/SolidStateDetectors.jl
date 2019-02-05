abstract type AbstractGeometry{T} end

function in(pt::AnyPoint{T}, vg::Vector{AbstractGeometry{T}})::Bool where {T <: AbstractFloat}
    is_inside::Bool = false
    for g in vg
        if in(pt, g) 
            is_inside = true
            break
        end
    end
    return is_inside
end

abstract type Volume{T} <: AbstractGeometry{T} end


struct Tubs{T <: AbstractFloat} <: Volume{T}
    name::String
    hierarchy::Int
    ϵ
    behaviour::String
    potential
    rStart::T
    rStop::T
    θStart::T
    θStop::T
    zStart::T
    zStop::T
    function Tubs(name::String, hierarchy::Int, ϵ, behaviour::String, potential, rStart::T, rStop::T, θStart::T, θStop::T, zStart::T, zStop::T) where T<:AbstractFloat
    return new{T}(name, hierarchy, ϵ, behaviour, potential, rStart, rStop, θStart, θStop, zStart, zStop)
    end
end

function Tubs(name, hierarchy, ϵ, behaviour, potential, rStart::T, rStop::T, θStart::T, θStop::T, zStart::T, zStop::T) where T<:AbstractFloat
    return Tubs{typeof(zStop)}(name, hierarchy, ϵ, behaviour, potential, rStart, rStop, θStart, θStop, zStart, zStop)
end


mutable struct CartesianBox3D{T} <: AbstractGeometry{T}
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

function Geometry(T::DataType, geometries::Vector, inputunit::Unitful.Units)#::AbstractGeometry{T} where {T <: AbstractFloat}
    return [Geometry(T, Val{Symbol(geometry["type"])}(), geometry, inputunit) for geometry in geometries]
end