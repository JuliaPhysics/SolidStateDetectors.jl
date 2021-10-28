function vertices(p::AbstractVolumePrimitive{T}, spacing::T)::Vector{CartesianPoint{T}} where {T}
    fs = surfaces(p)
    vs = Vector{CartesianPoint{T}}()
    for s in fs
        append!(vs, vertices(s, spacing))
    end
    vs
end

extremum(p::AbstractVolumePrimitive) = maximum([extremum(s) + norm(s.origin) for s in surfaces(p)])

@recipe function f(p::AbstractVolumePrimitive; n_samples = 200)
    fs = surfaces(p)
    # In principle it would be better to directly get the edges of the primitive
    # Can be improved later...
    linecolor --> :black
    if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        aspect_ratio --> 1.0
    end 
    @series begin
        label --> "$(nameof(typeof(p)))"
        if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :samplesurface
            vertices(p, extremum(p)/n_samples)
        else
            [fs...]
        end
    end
end