function vertices(p::AbstractVolumePrimitive{T}, spacing::T)::Vector{CartesianPoint{T}} where {T}
    fs = surfaces(p)
    vs = Vector{CartesianPoint{T}}()
    for s in fs
        append!(vs, vertices(s, spacing))
    end
    vs
end

@recipe function f(p::AbstractVolumePrimitive)
    fs = surfaces(p)
    # In principle it would be better to directly get the edges of the primitive
    # Can be improved later...
    linecolor --> :black
    if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        aspect_ratio --> 1.0
    end 
    @series begin
        label --> "$(nameof(typeof(p)))"
        [fs...]
    end
end

