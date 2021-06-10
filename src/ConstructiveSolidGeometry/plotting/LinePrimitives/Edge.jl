@recipe function f(e::Edge)
    linecolor --> :black
    @series begin
        seriestype --> :path3d
        [e.a[1], e.b[1]], [e.a[2], e.b[2]], [e.a[3], e.b[3]]
    end
end


