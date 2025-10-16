@recipe function f(pt::CartesianPoint)
    xguide --> "x"
    yguide --> "y"
    zguide --> "z"
    if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        aspect_ratio --> 1.0
    end
    @series begin
        seriestype --> :scatter
        [pt.x], [pt.y], [pt.z]
    end
end

@recipe function f(pt::CylindricalPoint)
    @series begin
        CartesianPoint(pt)
    end
end

@recipe function f(v::AbstractVector{<:CartesianPoint})
    default_markersize = 4.0
    n_thresh = 10
    n = length(v)
    invmarkerscale = sqrt(max(1, 2 + (n - 1 - n_thresh) / n_thresh))
    auto_markersize = max(default_markersize / invmarkerscale, 0.1)
    auto_markerstrokewidth = auto_markersize < 1 ? zero(auto_markersize) : 0.1 * default_markersize

    xguide --> "x"
    yguide --> "y"
    zguide --> "z"
    markersize --> auto_markersize
    markerstrokewidth --> auto_markerstrokewidth
    if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        aspect_ratio --> 1.0
    end 
    @series begin
        if !isempty(v) seriestype --> :scatter end
        [v[i].x for i in eachindex(v)], [v[i].y for i in eachindex(v)], [v[i].z for i in eachindex(v)]
    end
end

@recipe function f(v::AbstractVector{<:CylindricalPoint})
    @series begin
        map(x -> CartesianPoint(x), v)
    end
end
