function edges(e::Ellipse{T,T,Nothing}; n = 4) where {T}
    φs = range(0, stop = 2π, length = n + 1)
    pts = [CartesianPoint(CylindricalPoint{T}(e.r, φ, zero(T))) for φ in φs]
    pts = map(p -> _transform_into_global_coordinate_system(p, e), pts)
    edges = [Edge(pts[i], pts[i+1]) for i in 1:n]
end
function edges(e::Ellipse{T,T,Tuple{T,T}}; n = 4) where {T}
    φs = range(e.φ[1], stop = e.φ[2], length = n + 1)
    pts = [CartesianPoint(CylindricalPoint{T}(e.r, φ, zero(T))) for φ in φs]
    pts = map(p -> _transform_into_global_coordinate_system(p, e), pts)
    edges = [Edge(pts[i], pts[i+1]) for i in 1:n]
end
function edges(e::Ellipse{T,NTuple{2,Tuple{T}},Nothing}; n = 4) where {T}
    φs = range(0, stop = 2π, length = n + 1)
    pts = [CartesianPoint{T}(e.r[1][1] * cos(φ), e.r[2][1] * sin(φ), zero(T)) for φ in φs]
    pts = map(p -> _transform_into_global_coordinate_system(p, e), pts)
    edges = [Edge(pts[i], pts[i+1]) for i in 1:n]
end

@recipe function f(e::Ellipse; n = 40)
    @series begin
        label --> "Ellipse"
        edges(e, n = n)
    end
end


