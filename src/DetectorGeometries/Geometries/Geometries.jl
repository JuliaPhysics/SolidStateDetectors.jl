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

abstract type Volume{T} <: AbstractGeometry{T, 3} end


struct Tubs{T <: AbstractFloat} <: Volume{T}
    name::String
    hierarchy::Int
    ϵ
    behaviour::String
    potential
    rStart::T
    rStop::T
    φStart::T
    φStop::T
    zStart::T
    zStop::T
    function Tubs(name::String, hierarchy::Int, ϵ, behaviour::String, potential, rStart::T, rStop::T, φStart::T, φStop::T, zStart::T, zStop::T) where T<:AbstractFloat
    return new{T}(name, hierarchy, ϵ, behaviour, potential, rStart, rStop, φStart, φStop, zStart, zStop)
    end
end

function Tubs(name, hierarchy, ϵ, behaviour, potential, rStart::T, rStop::T, φStart::T, φStop::T, zStart::T, zStop::T) where T<:AbstractFloat
    return Tubs{typeof(zStop)}(name, hierarchy, ϵ, behaviour, potential, rStart, rStop, φStart, φStop, zStart, zStop)
end


