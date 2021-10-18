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
    edges::Vector{Vector{Int64}}
end

function mesh(p::AbstractSurfacePrimitive{T}; n_arc = 40, n_vert_lines = 2)::Mesh{T} where {T}
    vs = vertices(p, n_arc)
    x, y, z = broadcast(i -> getindex.(vs, i), (1,2,3))
    c = connections(p, n_arc)
    e = mesh_edges(p, n_arc, n_vert_lines)
    Mesh{T}(x,y,z,c,e)
end

@recipe function f(m::Mesh{T}; edgewidth = 0) where {T}
    seriestype --> :mesh3d
    xguide --> "x"
    yguide --> "y"
    zguide --> "z"
    unitformat --> :slash
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :mesh3d   
        @series begin
            if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
                linewidth --> 0.75
            else
                linewidth --> 0.1
            end
            linecolor --> :white
            seriescolor --> 1
            seriesalpha --> 0.5
            connections := m.connections
            m.x*internal_length_unit, m.y*internal_length_unit, m.z*internal_length_unit
        end
    end
    if edgewidth != 0 || (haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :path)
        @series begin
            label := ""
            seriescolor --> :black
            if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :path
                linewidth --> 1
            else
                linewidth := edgewidth
            end
            seriestype := :path
            [m.x[e]*internal_length_unit for e in m.edges], [m.y[e]*internal_length_unit for e in m.edges], [m.z[e]*internal_length_unit for e in m.edges]
        end
    end
end