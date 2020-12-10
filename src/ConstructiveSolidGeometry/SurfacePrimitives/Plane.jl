struct Plane{T,TP,TA} <: AbstractSurfacePrimitive{T}
    p1::TP
    p2::TP
    p3::TP
    add_points::TA
    function Plane( ::Type{T},
                   p1::CartesianPoint{T},
                   p2::CartesianPoint{T},
                   p3::CartesianPoint{T},
                   add_points:: Union{Nothing, CartesianPoint{T}, Vector{CartesianPoint{T}}}) where {T}
        new{T, CartesianPoint{T}, typeof(add_points)}(p1, p2, p3, add_points)
    end
end

#Constructors
function Plane(; p1 = CartesianPoint{Float32}(0,0,0), p2 = CartesianPoint{Float32}(1,0,0), p3 = CartesianPoint{Float32}(0,1,0), add_points = nothing)
    T = T = float(promote_type(typeof.((p1.x,p2.x,p3.x))...))
    Plane(T, p1, p2, p3, add_points)
end

Plane(p1, p2, p3, add_points) = Plane(;p1, p2, p3, add_points)

#plotting
function get_plot_points(p::Plane{T, CartesianPoint{T}, Nothing}; n = 30) where {T <: AbstractFloat}
    plot_points = Vector{CartesianPoint{T}}[]
    points = [p.p1, p.p2, p.p3]
    for i in 1:length(points)
        push!(plot_points, Vector{CartesianPoint{T}}([points[i], points[i%length(points)+1]]))
    end
    plot_points
end

function mesh(p::Plane{T}; n = 30) where {T <: AbstractFloat}
    points = [p.p1, p.p2, p.p3]
    Mesh(map(p -> p.x, points),map(p -> p.y, points),map(p -> p.z, points))
end
