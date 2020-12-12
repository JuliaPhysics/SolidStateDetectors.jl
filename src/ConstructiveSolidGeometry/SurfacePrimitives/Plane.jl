struct Plane{T,TP,TA} <: AbstractSurfacePrimitive{T}
    p1::TP
    p2::TP
    p3::TP
    p4::TA
    function Plane( ::Type{T},
                   p1::CartesianPoint{T},
                   p2::CartesianPoint{T},
                   p3::CartesianPoint{T},
                   p4:: Union{Nothing, CartesianPoint{T}}) where {T}
        new{T, CartesianPoint{T}, typeof(p4)}(p1, p2, p3, p4)
    end
end

#Constructors
function Triangle(;p1 = CartesianPoint{Float32}(0,0,0), p2 = CartesianPoint{Float32}(1,0,0), p3 = CartesianPoint{Float32}(0,1,0))
    T = float(promote_type(eltype.((p1, p2, p3))...))
    Plane(T, p1, p2, p3, nothing)
end

Triangle(p1, p2, p3) = Triangle(;p1, p2, p3)

function surface_vector(plane::Plane{T})::CartesianVector{T} where {T}
    v1::CartesianVector{T} = CartesianVector{T}(plane.p2 - plane.p1)
    v2::CartesianVector{T} = CartesianVector{T}(plane.p3 - plane.p1)
    normalize(cross(v1,v2))
end

function on_infinite_plane(p::CartesianPoint{T}, plane::Plane{T})::Bool where {T}
    v::CartesianVector{T} = surface_vector(plane)
    v1 = CartesianVector{T}(plane.p2 - p)
    v2 = CartesianVector{T}(plane.p3 - p)
    cross(v, cross(v1,v2)) ≈ CartesianVector{T}(0,0,0)
end

#returns points on infinite plane to be added to the plane
get_points_on_infinite_plane!(p::Nothing, plane::Plane{T}) where {T} = nothing

get_points_on_infinite_plane!(p::CartesianPoint{T}, plane::Plane{T}) where {T} = on_infinite_plane(p, plane) ? p : nothing

function get_points_on_infinite_plane!(vp::Vector{CartesianPoint{T}}, plane::Plane{T}) where {T}
    vp_onplane = CartesianPoint{T}[]
    for p in vp
        on_infinite_plane(p, plane) ? push!(vp_onplane, p) : nothing
    end
    length(vp_onplane) == 0 ? nothing : vp_onplane
end

function Plane(;p1 = CartesianPoint{Float32}(0,0,0), p2 = CartesianPoint{Float32}(1,0,0), p3 = CartesianPoint{Float32}(0,1,0), p4 = nothing)
    T = float(promote_type(eltype.((p1, p2, p3))...))
    triangle = Triangle(p1, p2, p3)
    Plane(T, p1, p2, p3, get_points_on_infinite_plane!(p4, triangle))
end

Plane(p1, p2, p3, p4) = Plane(;p1, p2, p3, p4)
Plane(p1, p2, p3) = Triangle(;p1, p2, p3)

get_vertices(plane::Plane{T, CartesianPoint{T}, Nothing}) where {T} = [plane.p1, plane.p2, plane.p3]
get_vertices(plane::Plane{T, CartesianPoint{T}, CartesianPoint{T}}) where {T} = [plane.p1, plane.p2, plane.p3, plane.p4]


function sample(plane::Plane{T, CartesianPoint{T}, Nothing}, step::R) where {T, R<:Real}
    sampled_points = CartesianPoint{T}[]
    v1::CartesianVector{T} = CartesianVector{T}(plane.p2 - plane.p1)
    v2::CartesianVector{T} = CartesianVector{T}(plane.p3 - plane.p1)
    #u = range(0,1, length = n)
    #v = range(0,1, length = n)
    step_u = step/norm(v1)
    step_v = step/norm(v2)
    u = 0:step_u:1
    v = 0:step_v:1
    for u_i in u
        j = 1
        while j <= length(v) && v[j] + u_i <= 1
            v_t = CartesianPoint{T}(u_i*v1 + v[j]*v2)
            push!(sampled_points, v_t)
            j = j + 1
        end
    end
    sampled_points
end

#plotting

function get_plot_points(plane::Plane{T}; n = 30) where {T <: AbstractFloat}
    plot_points = Vector{CartesianPoint{T}}[]
    vertices = get_vertices(plane)
    for i in 1:length(vertices)
        push!(plot_points, Vector{CartesianPoint{T}}([vertices[i], vertices[i%length(vertices)+1]]))
    end
    plot_points
end

function mesh(plane::Plane{T, CartesianPoint{T}, Nothing}; n = 30) where {T <: AbstractFloat}
    vertices = get_vertices(plane)
    Mesh(map(p -> p.x, vertices),map(p -> p.y, vertices),map(p -> p.z, vertices))
end
