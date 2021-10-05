@recipe function f(pt::CartesianPoint)
    xguide --> "x"
    yguide --> "y"
    zguide --> "z"
    unitformat --> :slash
    if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        aspect_ratio --> 1.0
    end 
    @series begin
        seriesstyle --> :scatter
        [pt.x]u"m", [pt.y]u"m", [pt.z]u"m"
    end
end

@recipe function f(pt::CylindricalPoint)
    @series begin
        CartesianPoint(pt)
    end
end

@recipe function f(v::AbstractVector{<:CartesianPoint})
    xguide --> "x"
    yguide --> "y"
    zguide --> "z"
    unitformat --> :slash
    if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        aspect_ratio --> 1.0
    end 
    @series begin
        seriesstyle --> :scatter
        [v[i].x for i in eachindex(v)]u"m", [v[i].y for i in eachindex(v)]u"m", [v[i].z for i in eachindex(v)]u"m"
    end
end

@recipe function f(v::AbstractVector{<:CylindricalPoint})
    @series begin
        map(x -> CartesianPoint(x), v)
    end
end