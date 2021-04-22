function get_plot_points(c::ConalPlane{T}; n = 30) where {T <: AbstractFloat}
    v = get_vertices(c)
    [Vector{CartesianPoint{T}}([v[1], v[2]]), Vector{CartesianPoint{T}}([v[2], v[3]]), Vector{CartesianPoint{T}}([v[3], v[4]]), Vector{CartesianPoint{T}}([v[4], v[1]])]
end

function mesh(c::ConalPlane{T}; n = 30) where {T <: AbstractFloat}
    v = unique(get_vertices(c))
    while length(v) < 3
        push!(v,v[1])
    end
    mesh(Plane(v..., p4_on_plane_check = false))
 end
