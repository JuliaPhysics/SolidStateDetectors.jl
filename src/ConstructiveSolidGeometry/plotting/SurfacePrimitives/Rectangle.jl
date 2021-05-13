function get_plot_points(r::Rectangle{T}; n = 30) where {T <: AbstractFloat}
    v = get_vertices(r)
    [Vector{CartesianPoint{T}}([v[1], v[2]]), Vector{CartesianPoint{T}}([v[2], v[4]]), Vector{CartesianPoint{T}}([v[3], v[4]]), Vector{CartesianPoint{T}}([v[3], v[1]])]
end

function mesh(r::Rectangle{T}; n = 30) where {T <: AbstractFloat}
    vertices = [get_vertices(r)...]
    Mesh(reshape(map(p -> p.x, vertices), (2,2)),reshape(map(p -> p.y, vertices), (2,2)),reshape(map(p -> p.z, vertices), (2,2)))
 end
