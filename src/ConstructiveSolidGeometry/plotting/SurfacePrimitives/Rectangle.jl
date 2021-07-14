function get_plot_points(r::Rectangle{<:Any, T}; n = 30) where {T <: AbstractFloat}
    v = get_vertices(r)
    [Vector{CartesianPoint{T}}([v[1], v[2]]), Vector{CartesianPoint{T}}([v[2], v[4]]), Vector{CartesianPoint{T}}([v[3], v[4]]), Vector{CartesianPoint{T}}([v[3], v[1]])]
end

function mesh(r::Rectangle{<:Any, T}; n = 30) where {T <: AbstractFloat}
    vertices = [get_vertices(r)...]
    Mesh(reshape(broadcast(p -> p.x, vertices), (2,2)),reshape(broadcast(p -> p.y, vertices), (2,2)),reshape(broadcast(p -> p.z, vertices), (2,2)))
 end
