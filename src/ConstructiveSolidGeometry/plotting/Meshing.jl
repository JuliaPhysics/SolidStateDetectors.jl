"""
    struct Mesh{T}

* `x`: X-coordinate mesh in meters
* `y`: Y-coordinate mesh in meters
* `z`: Z-coordinate mesh in meters
* `connections`: defines which points are connected to eah other

"""

struct Mesh{T}
    x::Array{T}
    y::Array{T}
    z::Array{T}
    connections::Vector{Vector{Int64}} 
end

function rotate(mesh::Mesh{T}, R::AbstractMatrix)::Mesh{T} where {T}
    n = length(mesh.x)
    x = zeros(n)
    y = zeros(n)
    z = zeros(n)
    for i in 1:n
        v = [mesh.x[i], mesh.y[i], mesh.z[i]]
        v = R*v
        x[i], y[i], z[i] = v[1], v[2], v[3]
    end
    Mesh{T}(x,y,z,mesh.connections)
end
(*)(R::AbstractMatrix, mesh::Mesh{T}) where {T} = rotate(mesh, R)

translate(mesh::Mesh{T}, t::CartesianPoint) where {T} =
    Mesh{T}(mesh.x .+ t.x, mesh.y .+ t.y, mesh.z .+ t.z, mesh.connections)
(+)(mesh::Mesh, t::CartesianPoint) = translate(mesh,t)

@recipe function f(m::Mesh{T}) where {T}
    seriestype := :mesh3d
    linewidth --> 0.1
    linecolor --> :white
    seriescolor --> 1
    seriesalpha --> 0.5
    connections := m.connections
    m.x, m.y, m.z
end