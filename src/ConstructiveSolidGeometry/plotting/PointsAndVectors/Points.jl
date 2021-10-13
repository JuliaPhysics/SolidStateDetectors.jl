@recipe function f(pt::CartesianPoint)
    xguide --> "x"
    yguide --> "y"
    zguide --> "z"
    unitformat --> :slash
    if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        aspect_ratio --> 1.0
    end 
    @series begin
        seriestype --> :scatter
        [pt.x]*internal_length_unit, [pt.y]*internal_length_unit, [pt.z]*internal_length_unit
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
        seriestype --> :scatter
        [v[i].x for i in eachindex(v)]*internal_length_unit, [v[i].y for i in eachindex(v)]*internal_length_unit, [v[i].z for i in eachindex(v)]*internal_length_unit
    end
end

@recipe function f(v::AbstractVector{<:CylindricalPoint})
    @series begin
        map(x -> CartesianPoint(x), v)
    end
end