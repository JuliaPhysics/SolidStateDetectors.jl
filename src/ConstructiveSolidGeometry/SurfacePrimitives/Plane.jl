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

Triangle(p1, p2, p3) = Triangle(;p1 = p1, p2 = p2, p3 = p3)

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
    abs(dot(n,v))
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

function distance_to_line_segment(point::AbstractCoordinatePoint{T}, seg::Tuple{AbstractCoordinatePoint{T},AbstractCoordinatePoint{T}})::T where {T}
    point = CartesianPoint(point)
    seg = (CartesianPoint(seg[1]), CartesianPoint(seg[2]))
    v12 = normalize(CartesianVector{T}(seg[2] - seg[1]))
    v_point_1 = CartesianVector{T}(point - seg[1])
    proj_on_v12 = dot(v12,v_point_1)
    if geom_round(proj_on_v12) ≤ T(0)
        return norm(seg[1] - point)
    else
        v_point_2 = CartesianVector{T}(point - seg[2])
        if geom_round(dot(v12,v_point_2)) ≥ T(0)
            return norm(seg[2] - point)
        else
            return sqrt(abs(dot(v_point_1,v_point_1) - proj_on_v12^2))
        end
    end
end

function distance_to_surface(point::AbstractCoordinatePoint{T}, tri::Plane{T, CartesianPoint{T}, Nothing})::T where {T}
    point = CartesianPoint(point)
    u, v = get_planar_coordinates(point, tri)
    if geom_round(u) ≥ T(0) && geom_round(v) ≥ T(0) #++ quadrant. Origin is tri.p1
        if geom_round(u+v) ≤ T(1) #in triangle
            return distance_to_infinite_plane(point, tri)
        else # on side 2_3 of tri
            return distance_to_line_segment(point, (tri.p2,tri.p3))
        end
    elseif geom_round(u) ≤ T(0) && geom_round(v) ≤ T(0) #-- quadrant
        return norm(tri.p1 - point)
    elseif geom_round(u) > T(0) && geom_round(v) < T(0) #+- quadrant, on side 1_2 of tri
        return distance_to_line_segment(point, (tri.p1,tri.p2))
    elseif geom_round(u) < T(0) && geom_round(v) > T(0) #-+ quadrant, on side 3_1 of tri
        return distance_to_line_segment(point, (tri.p3,tri.p1))
    end
end

function sample(tri::Plane{T, CartesianPoint{T}, Nothing}, step::Real) where {T}
    v1, v2 = get_spanning_vectors(tri)
    step_u = step/norm(v1)
    step_v = step/norm(v2)
    samples = [
        CartesianPoint{T}(tri.p1 + u*v1 + v*v2)
        for u in 0:step_u:1
        for v in 0:step_v:1-u
    ]
end

function sample(tri::Plane{T, CartesianPoint{T}, Nothing}, Nsamps::NTuple{3,Int}) where {T}
    v1, v2 = get_spanning_vectors(tri)
    samples = [
        CartesianPoint{T}(tri.p1 + u*v1 + v*v2)
        for u in (Nsamps[1] ≤ 1 ? 0 : range(0, 1, length = Nsamps[1]))
        for v in (Nsamps[2] ≤ 1 ? 0 : range(0, 1-u, length = Nsamps[2]))
    ]
end

function Quadrilateral(;p1 = CartesianPoint{Float32}(0,0,0), p2 = CartesianPoint{Float32}(1,0,0), p3 = CartesianPoint{Float32}(1,1,0), p4 = CartesianPoint{Float32}(0,1,0), p4_on_plane_check = true)
    T = float(promote_type(eltype.((p1, p2, p3, p4))...))
    tri = Triangle(p1, p2, p3)
    #will return triangle if conditions are not met. Order of points matters, a continous non intersecting line needs to be drawn in p1->p2->p3->p4->p1
    if geom_round(p4) in [geom_round(p1), geom_round(p2), geom_round(p3)]
        println("Identical vertex")
        return tri
    elseif p4_on_plane_check
        if on_infinite_plane(p4, tri)
            if !projection_in_triangle(p4, tri)
                v = CartesianVector{T}(p4 - p1)
                v2 = CartesianVector{T}(p3 - p1)
                n = get_surface_vector_nonunitary(tri)
                if geom_round(dot(n,cross(v2,v))) ≥ 0
                    return Plane(T, p1, p2, p3, p4)
                else #intersecting segments
                    return tri
                end
            else
                return Plane(T, p1, p2, p3, p4)
            end
        else #not on plane
            return tri
        end
    else
        return Plane(T, p1, p2, p3, p4)
    end
end

Quadrilateral(p1, p2, p3, p4; p4_on_plane_check = true) = Quadrilateral(;p1 = p1, p2 = p2, p3 = p3, p4 = p4, p4_on_plane_check = p4_on_plane_check)
Plane(p1, p2, p3, p4; p4_on_plane_check = true) = Quadrilateral(;p1 = p1, p2 = p2, p3 = p3, p4 = p4, p4_on_plane_check = p4_on_plane_check)
Plane(p1, p2, p3; p4_on_plane_check = true) = Triangle(;p1 = p1, p2 = p2, p3 = p3)

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

function sample(quad::Plane{T, CartesianPoint{T}, CartesianPoint{T}}, sampling::Union{Real, NTuple{3,Int}}) where {T}
    tri1, tri2 = decompose_into_tiangles(quad)
    samples = sample(tri1, sampling)
    append!(samples, sample(tri2, sampling))
end
