"""
    struct Mesh{T}

* `x`: X-coordinate mesh in meters
* `y`: Y-coordinate mesh in meters
* `z`: Z-coordinate mesh in meters

(X[i,j], Y[i,j], Z[i,j]) is a cartesian point
"""

struct Mesh{T}
    x::Array{T,2}
    y::Array{T,2}
    z::Array{T,2}
    function Mesh(
        x::Array{T,2},
        y::Array{T,2},
        z::Array{T,2}) where {T}
        @assert size(x) == size(y) && size(x) == size(z)
        new{T}(x, y, z)
    end
end

function get_cartesian_point_from_mesh(m::Mesh{T}, index::Tuple{I,I}) where {T, I<:Int}
    i, j = index
    CartesianPoint{T}(m.x[i,j], m.y[i,j], m.z[i,j])
end

size(m::Mesh{T}) where {T} = size(m.x)

function get_plot_meshes(v::AbstractVolumePrimitive{T}) where {T <: AbstractFloat}
  surfaces = get_decomposed_surfaces(v)
  meshes = Mesh{T}[]
  for surf in surfaces
      push!(meshes, mesh(surf))
  end
  meshes
end

function get_plot_meshes(s::AbstractSurfacePrimitive{T}) where {T <: AbstractFloat}
  meshes = Mesh{T}[]
  push!(meshes, mesh(s))
  meshes
end
