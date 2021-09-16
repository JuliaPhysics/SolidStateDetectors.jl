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
end

struct PolyMesh{T,N}
    x::Vector{SVector{N,T}}
    y::Vector{SVector{N,T}}
    z::Vector{SVector{N,T}}
end

function get_cartesian_point_from_mesh(m::Mesh{T}, index::Tuple{I,I}) where {T, I<:Int}
    i, j = index
    CartesianPoint{T}(m.x[i,j], m.y[i,j], m.z[i,j])
end

size(m::Mesh{T}) where {T} = size(m.x)

function rotate(mesh::Mesh{T}, R::AbstractMatrix)::Mesh{T} where {T}
    n, m = size(mesh)
    X = zero(mesh.x)
    Y = zero(mesh.y)
    Z = zero(mesh.z)
    for i in 1:n
        for j in 1:m
            v = get_cartesian_point_from_mesh(mesh, (i,j))
            v = R*v
            X[i,j], Y[i,j], Z[i,j] = v[1], v[2], v[3]
        end
    end
    Mesh{T}(X,Y,Z)
end
(*)(R::AbstractMatrix, mesh::Mesh{T}) where {T} = rotate(mesh, R)

translate(mesh::Mesh{T}, t::CartesianPoint) where {T} =
    Mesh{T}(mesh.x .+ t.x, mesh.y .+ t.y, mesh.z .+ t.z)
(+)(mesh::Mesh, t::CartesianPoint) = translate(mesh,t)

function polymesh(mesh::Mesh{T})::PolyMesh{T,4} where {T}
    #mesh = RotXY(0.0000001,0.0000001) * mesh
    n, m = size(mesh)
    x = [reshape(mesh.x[i:i+1,j:j+1], 4) for i in 1:n-1, j in 1:m-1]
    y = [reshape(mesh.y[i:i+1,j:j+1], 4) for i in 1:n-1, j in 1:m-1]
    z = [reshape(mesh.z[i:i+1,j:j+1], 4) for i in 1:n-1, j in 1:m-1]
    dim = (n-1)*(m-1)
    PolyMesh{T,4}(reshape(x,dim), reshape(y,dim), reshape(z,dim))
end

#=
function gr_mesh_patch!(x::Vector{Vector{T}}, y::Vector{Vector{T}}, z::Vector{Vector{T}}) where {T}
    for i in 1:length(x)
        append!(x[i], sum(x[i])/4)
        append!(y[i], sum(y[i])/4)
        append!(z[i], sum(z[i])/4)
    end
end
=#

@recipe function f(m::PolyMesh{T,4}) where {T} #only supported in gr()
    seriestype := :mesh3d
    colorbar := false
    seriescolor --> 1
    seriesalpha --> 0.5
    connections := ([1,0],[0,2],[3,3])
    @series begin
        m.x[1], m.y[1], m.z[1]
    end
    for i in 2:length(m.x)
        @series begin
            label := nothing
            m.x[i], m.y[i], m.z[i]
        end
    end
end

@recipe function f(m::Mesh{T}) where {T}
    seriestype := :surface
    colorbar := false
    linewidth --> 0.1
    linecolor --> :white
    seriescolor --> [1,1]
    seriesalpha --> 0.5
    m.x, m.y, m.z
end