@recipe function f(e::Edge)
    linecolor --> :black
    label --> "Edge"
    xguide --> "x"
    yguide --> "y"
    zguide --> "z"
    unitformat --> :slash
    if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        aspect_ratio --> 1.0
    end 
    @series begin
        seriestype --> :path3d
        [e.a[1], e.b[1]]u"m", [e.a[2], e.b[2]]u"m", [e.a[3], e.b[3]]u"m"
    end
end


