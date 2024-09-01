function sample(p::AbstractVolumePrimitive{T}, spacing::T)::Vector{CartesianPoint{T}} where {T}
    fs = surfaces(p)
    vs = Vector{CartesianPoint{T}}()
    for s in fs
        append!(vs, sample(s, spacing))
    end
    vs
end

@recipe function f(p::AbstractVolumePrimitive{T}; n_samples = 40, slice_val = T(0)) where {T}
    linecolor --> :black
    seriestype --> :csg
    @series begin
        label --> "$(nameof(typeof(p)))"
        projections = [:x, :y, :z]
        if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :samplesurface
            sample(p, extremum(p)/n_samples)
        elseif haskey(plotattributes, :seriestype) && plotattributes[:seriestype] in projections
            #Unlike f(CartesianPoint), there is no ssd plot recipe for f(array,array) so must define attributes here
            xguide --> "x"
            xunit --> internal_length_unit
            yguide --> "y"
            yunit --> internal_length_unit
            unitformat --> :slash
            spacing = extremum(p)/n_samples
            samples = filter(pt -> abs(getproperty(pt, plotattributes[:seriestype]) - slice_val) < spacing/2, sample(p, spacing))
            proj = filter(x -> x != plotattributes[:seriestype], projections)
            internal_length_unit*getproperty.(samples, proj[1]), internal_length_unit*getproperty.(samples, proj[2])
        else
            [surfaces(p)...]
        end
    end
end