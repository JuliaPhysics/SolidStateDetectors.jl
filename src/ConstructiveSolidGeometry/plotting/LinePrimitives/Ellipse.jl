function edges(e::Ellipse{T,T,Nothing}; n = 4) where {T}
    φs = range(0, stop = 2π, length = n + 1)
    pts = [CartesianPoint(CylindricalPoint{T}(e.r, φ, zero(T))) for φ in φs]
    pts = map(p -> _transform_into_global_coordinate_system(p, e), pts)
    edges = [Edge(pts[i], pts[i+1]) for i in 1:n]
end

@recipe function f(e::Ellipse; n = 40)
    @series begin
        label --> "Ellipse"
        edges(e, n = n)
    end
end


