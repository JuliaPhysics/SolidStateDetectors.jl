@recipe function f(pt::CartesianPoint)
    xguide --> "X"
    yguide --> "Y"
    zguide --> "Z"
    @series begin
        seriesstyle --> :scatter
        [pt.x], [pt.y], [pt.z]
    end
end

@recipe function f(pt::CylindricalPoint)
    @series begin
        CartesianPoint(pt)
    end
end

@recipe function f(v::AbstractVector{<:CartesianPoint})
    xguide --> "X"
    yguide --> "Y"
    zguide --> "Z"
    @series begin
        seriesstyle --> :scatter
        [v[i].x for i in eachindex(v)], [v[i].y for i in eachindex(v)], [v[i].z for i in eachindex(v)]
    end
end

@recipe function f(v::AbstractVector{<:CylindricalPoint})
    @series begin
        map(x -> CartesianPoint(x), v)
    end
end