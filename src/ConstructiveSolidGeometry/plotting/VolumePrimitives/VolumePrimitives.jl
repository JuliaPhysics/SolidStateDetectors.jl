function sample(p::AbstractVolumePrimitive{T}, spacing::T)::Vector{CartesianPoint{T}} where {T}
    fs = surfaces(p)
    vs = Vector{CartesianPoint{T}}()
    for s in fs
        append!(vs, sample(s, spacing))
    end
    vs
end

@recipe function f(p::AbstractVolumePrimitive; n_samples = 40)
    linecolor --> :black
    if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        aspect_ratio --> 1.0
    end 
    @series begin
        label --> "$(nameof(typeof(p)))"
        if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :samplesurface
            sample(p, extremum(p)/n_samples)
        else
            [surfaces(p)...]
        end
    end
end