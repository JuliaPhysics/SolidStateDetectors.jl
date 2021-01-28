function get_plot_points(plane::Plane{T}; n = 30) where {T <: AbstractFloat}
    plot_points = Vector{CartesianPoint{T}}[]
    vertices = get_vertices(plane)
    for i in 1:length(vertices)
        push!(plot_points, Vector{CartesianPoint{T}}([vertices[i], vertices[i%length(vertices)+1]]))
    end
    plot_points
end

function mesh(tri::Plane{T, CartesianPoint{T}, Nothing}; n = 30) where {T <: AbstractFloat}
    vertices = [tri.p1]
    append!(vertices, get_vertices(tri))
    Mesh(reshape(map(p -> p.x, vertices), (2,2)),reshape(map(p -> p.y, vertices), (2,2)),reshape(map(p -> p.z, vertices), (2,2)))
end

function mesh(quad::Plane{T, CartesianPoint{T}, CartesianPoint{T}}; n = 30) where {T <: AbstractFloat}
    tri1, tri2 = decompose_into_tiangles(quad)
    vertices = get_vertices(tri1)
    push!(vertices, tri2.p1)
    Mesh(reshape(map(p -> p.x, vertices), (2,2)),reshape(map(p -> p.y, vertices), (2,2)),reshape(map(p -> p.z, vertices), (2,2)))
end
