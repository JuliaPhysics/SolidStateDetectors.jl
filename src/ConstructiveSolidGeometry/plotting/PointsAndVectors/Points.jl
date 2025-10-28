@recipe function f(pt::CartesianPoint)
    xguide --> "x"
    yguide --> "y"
    zguide --> "z"
    if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        aspect_ratio --> 1.0
    end
    @series begin
        seriestype --> :scatter
        [to_internal_units(pt.x) * internal_length_unit], 
        [to_internal_units(pt.y) * internal_length_unit], 
        [to_internal_units(pt.z) * internal_length_unit]
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
        [to_internal_units(v[i].x) * internal_length_unit for i in eachindex(v)],
        [to_internal_units(v[i].y) * internal_length_unit for i in eachindex(v)],
        [to_internal_units(v[i].z) * internal_length_unit for i in eachindex(v)]
    end
end

@recipe function f(v::AbstractVector{<:CylindricalPoint})
    @series begin
        map(x -> CartesianPoint(x), v)
    end
end
