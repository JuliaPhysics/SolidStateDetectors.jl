"""
    struct Mesh{T}

* `x`: X-coordinate mesh in meters
* `y`: Y-coordinate mesh in meters
* `z`: Z-coordinate mesh in meters
* `connections`: defines which points are connected to eah other

"""

struct Mesh{T}
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    connections::Vector{Vector{Int64}} 
end

function mesh(p::AbstractSurfacePrimitive{T}; n_arc = 40)::Mesh{T} where {T}
    vs = vertices(p, n_arc)
    x, y, z = broadcast(i -> getindex.(vs, i), (1,2,3))
    c = connections(p, n_arc)
    Mesh{T}(x,y,z,c)
end


@recipe function f(m::Mesh{T}) where {T}
    seriestype := :mesh3d
    linewidth --> 0.75
    linecolor --> :white
    seriescolor --> 1
    seriesalpha --> 0.5
    connections := m.connections
    m.x, m.y, m.z
end