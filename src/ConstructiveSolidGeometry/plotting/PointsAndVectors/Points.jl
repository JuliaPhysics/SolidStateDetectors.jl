@recipe function f(p::CartesianPoint)
    @series begin
        seriesstyle --> :scatter
        [p.x], [p.y], [p.z]
    end
end

@recipe function f(p::CylindricalPoint)
    @series begin
        CartesianPoint(p)
    end
end

@recipe function f(v::AbstractVector{<:CartesianPoint})
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