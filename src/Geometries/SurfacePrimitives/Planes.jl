# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

# struct Plane{T} <: AbstractSurfacePrimitive{T, 3}
#     org::CartesianPoint{T}
#     normal::CartesianVector{T}
# end

# function in(pt::CartesianPoint{T}, p::Plane{T}; atol::Real = T(0))::Bool where {T}
#     if isapprox(pt, f.points[1], atol=atol) return true end
#     d::CartesianVector{T} = pt - f.points[1]
#     return -atol <= d ⋅ f.normal <= atol
# end

# function intersection(l::AbstractLine{T, 3, :cartesian}, f::Plane3D{T}; atol::Real = T(0)) where {T}
#     denom::T = f.normal ⋅ l.dir
#     if isapprox(denom, 0, atol = atol)
#         return false, CartesianPoint{T}(0, 0, 0)
#     else
#         u::T = f.normal ⋅ (f.points[1] - l.org) / denom
#         int_p::CartesianPoint{T} = l.org + u * l.dir 
#         return true, int_p
#     end
# end

# ToDo: use AbstractCoordinatePoint/Vector....
struct Plane
    reference_point::AbstractArray
    n⃗::AbstractArray

    function Plane(n⃗, reference_point)
        return new(n⃗, reference_point)
    end
end

function project_to_plane(v⃗::AbstractArray, n⃗::AbstractArray) #Vector to be projected, #normal vector of plane
    # solve (v⃗+λ*n⃗) ⋅ n⃗ = 0
    # n⃗ = n⃗ ./ norm(n⃗)
    λ = -1 * dot(v⃗, n⃗) / dot(n⃗, n⃗)
    SVector{3,eltype(v⃗)}(v⃗[1] + λ * n⃗[1], v⃗[2] + λ * n⃗[2], v⃗[3] + λ * n⃗[3])
end
