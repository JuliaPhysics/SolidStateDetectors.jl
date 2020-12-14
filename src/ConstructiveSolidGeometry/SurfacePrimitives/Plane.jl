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

get_vertices(tri::Plane{T, CartesianPoint{T}, Nothing}) where {T} = [tri.p1, tri.p2, tri.p3]

get_spanning_vectors(plane::Plane{T}) where {T} = (CartesianVector{T}(plane.p2 - plane.p1), CartesianVector{T}(plane.p3 - plane.p1))

function get_surface_vector(plane::Plane{T})::CartesianVector{T} where {T}
    normalize(cross(get_spanning_vectors(plane)...))
end

function get_surface_vector_nonunitary(plane::Plane{T})::CartesianVector{T} where {T}
    cross(get_spanning_vectors(plane)...)
end

function distance_to_infinite_plane(point::CartesianPoint{T}, plane::Plane{T})::T where {T}
    n::CartesianVector{T} = get_surface_vector(plane)
    v = CartesianVector{T}(point - plane.p1)
    geom_round(abs(dot(n,v)))
end

function on_infinite_plane(point::CartesianPoint{T}, plane::Plane{T})::Bool where {T}
    distance_to_infinite_plane(point, plane) == T(0)
end

function get_planar_coordinates(point::CartesianPoint{T}, tri::Plane{T, CartesianPoint{T}, Nothing})::Tuple{T,T} where {T}
    v1, v2 = get_spanning_vectors(tri)
    v = CartesianVector{T}(point - tri.p1)
    #equation to solve Ax = v where x is new 2D representation of plane. We solve by hitting it with A^T from right
    A = hcat(v1,v2)
    A_T = transpose(A)
    x = inv(A_T*A)*(A_T*v)
    x[1], x[2]
end

function projection_in_triangle(point::CartesianPoint{T}, tri::Plane{T, CartesianPoint{T}, Nothing})::Bool where {T}
    u, v = get_planar_coordinates(point, tri)
    geom_round(u) ≥ T(0) && geom_round(v) ≥ T(0) && geom_round(u+v) ≤ T(1)
end

in(point::CartesianPoint{T}, tri::Plane{T, CartesianPoint{T}, Nothing}) where {T} = projection_in_triangle(point, tri) && on_infinite_plane(point, tri)
#Order of checks makes no difference in performance in all tested cases

function _distance_to_line_segment(point::CartesianPoint{T}, seg::Tuple{CartesianPoint{T},CartesianPoint{T}})::T where {T}
    v12 = normalize(CartesianVector{T}(seg[2] - seg[1]))
    v_point_1 = CartesianVector{T}(point - seg[1])
    proj_on_v12 = dot(v12,v_point_1)
    if geom_round(proj_on_v12) ≤ T(0) #projection of point not on side 2->3
        return norm(seg[1] - point)
    else
        v_point_2 = CartesianVector{T}(point - seg[2])
        if geom_round(dot(v12,v_point_2)) ≥ T(0) #projection of point not on side 2->3
            return norm(seg[2] - point)
        else #projection of point on side 2->3
            return sqrt(dot(v_point_1,v_point_1) - proj_on_v12^2)
        end
    end
end

function distance_to_surface(point::CartesianPoint{T}, tri::Plane{T, CartesianPoint{T}, Nothing})::T where {T}
    u, v = get_planar_coordinates(point, tri)
    if geom_round(u) ≥ T(0) && geom_round(v) ≥ T(0) #++ quadrant. Origin is tri.p1
        if geom_round(u+v) ≤ T(1) #in triangle
            return distance_to_infinite_plane(point, tri)
        else # on side 2_3 of tri
            return _distance_to_line_segment(point, (tri.p2,tri.p3))
        end
    elseif geom_round(u) ≤ T(0) && geom_round(v) ≤ T(0) #-- quadrant
        return norm(tri.p1 - point)
    elseif geom_round(u) > T(0) && geom_round(v) < T(0) #+- quadrant, on side 1_2 of tri
        return _distance_to_line_segment(point, (tri.p1,tri.p2))
    elseif geom_round(u) < T(0) && geom_round(v) > T(0) #-+ quadrant, on side 3_1 of tri
        return _distance_to_line_segment(point, (tri.p3,tri.p1))
    end
end

function sample(tri::Plane{T, CartesianPoint{T}, Nothing}, step::Quantity{<:Real, Unitful.𝐋}) where {T}
    samples = CartesianPoint{T}[]
    v1, v2 = get_spanning_vectors(tri)
    #u = range(0,1, length = n)
    #v = range(0,1, length = n)
    step = T(ustrip(uconvert(u"m", step)))
    step_u = step/norm(v1)
    step_v = step/norm(v2)
    u = 0:step_u:1
    v = 0:step_v:1
    for u_i in u
        j = 1
        while j <= length(v) && v[j] + u_i <= 1
            v_t = CartesianPoint{T}(tri.p1 + u_i*v1 + v[j]*v2)
            push!(samples, v_t)
            j = j + 1
        end
    end
    samples
end

function Quadrilateral(;p1 = CartesianPoint{Float32}(0,0,0), p2 = CartesianPoint{Float32}(1,0,0), p3 = CartesianPoint{Float32}(1,1,0), p4 = CartesianPoint{Float32}(0,1,0))
    T = float(promote_type(eltype.((p1, p2, p3))...))
    tri = Triangle(p1, p2, p3)
    #will return triangle if conditions are not met
    if geom_round(p4) in [geom_round(p1), geom_round(p2), geom_round(p3)]
        return tri
    else
        if on_infinite_plane(p4, tri)
            if !projection_in_triangle(p4, tri)
                v = CartesianVector{T}(p4 - p1)
                v2 = CartesianVector{T}(p3 - p1)
                n = get_surface_vector_nonunitary(tri)
                if geom_round(dot(n,cross(v2,v))) ≥ 0
                    return Plane(T, p1, p2, p3, p4)
                else
                    return tri
                end
            else
                return Plane(T, p1, p2, p3, p4)
            end
        else
            return tri
        end
    end
end

Quadrilateral(p1, p2, p3, p4) = Quadrilateral(;p1, p2, p3, p4)
Plane(p1, p2, p3, p4) = Quadrilateral(;p1, p2, p3, p4)
Plane(p1, p2, p3) = Triangle(;p1, p2, p3)

get_vertices(quad::Plane{T, CartesianPoint{T}, CartesianPoint{T}}) where {T} = [quad.p1, quad.p2, quad.p3, quad.p4]


function decompose_into_tiangles(quad::Plane{T, CartesianPoint{T}, CartesianPoint{T}})::Tuple{Plane{T, CartesianPoint{T}, Nothing},Plane{T, CartesianPoint{T}, Nothing}} where {T}
    quad.p4 in Triangle(quad.p2,quad.p1,quad.p3) || quad.p2 in Triangle(quad.p4,quad.p1,quad.p3) ? (Triangle(quad.p1,quad.p2,quad.p4), Triangle(quad.p3,quad.p2,quad.p4)) : (Triangle(quad.p2,quad.p1,quad.p3), Triangle(quad.p4,quad.p1,quad.p3))
end


function in(point::CartesianPoint{T}, quad::Plane{T, CartesianPoint{T}, CartesianPoint{T}})::Bool where {T}
    tri1, tri2 = decompose_into_tiangles(quad)
    point in tri1 || point in tri2
end

function distance_to_surface(point::CartesianPoint{T}, quad::Plane{T, CartesianPoint{T}, CartesianPoint{T}})::T where {T}
    tri1, tri2 = decompose_into_tiangles(quad)
    min(distance_to_surface(point, tri1), distance_to_surface(point, tri2))
end

function sample(quad::Plane{T, CartesianPoint{T}, CartesianPoint{T}}, step::Quantity{<:Real, Unitful.𝐋}) where {T}
    tri1, tri2 = decompose_into_tiangles(quad)
    samples = sample(tri1, step)
    append!(samples, sample(tri2, step))
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
