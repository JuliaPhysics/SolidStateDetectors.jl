function get_plot_points(c::ConalPlane{T}; n = 30) where {T <: AbstractFloat}
    v = get_vertices(c)
    [Vector{CartesianPoint{T}}([v[1], v[2]]), Vector{CartesianPoint{T}}([v[2], v[4]]), Vector{CartesianPoint{T}}([v[3], v[4]]), Vector{CartesianPoint{T}}([v[3], v[1]])]
end

function mesh(c::ConalPlane{T}; n = 30) where {T <: AbstractFloat}
    vertices = unique(get_vertices(c))
    while length(vertices) < 4
        insert!(vertices, 1, vertices[1])
    end
    Mesh(reshape(map(p -> p.x, vertices), (2,2)),reshape(map(p -> p.y, vertices), (2,2)),reshape(map(p -> p.z, vertices), (2,2)))
 end
