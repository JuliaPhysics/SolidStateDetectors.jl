@recipe function f(t::Torus{T}; n = 30, seriescolor = :green) where {T}
    linewidth --> 2
    n --> n
    @series begin
        seriescolor --> seriescolor
        label --> "Torus"
        []
    end
    label := ""
    seriescolor := seriescolor
    LineSegments(t)
end
