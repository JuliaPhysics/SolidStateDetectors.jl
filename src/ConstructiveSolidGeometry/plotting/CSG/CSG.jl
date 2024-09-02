@inline primitives(vp::AbstractVolumePrimitive) = (vp,)
@inline function primitives(csg::AbstractConstructiveGeometry)
    (
        primitives(csg.a)...,
        primitives(csg.b)...
    )
end

function _get_csg_plot_objects(csg::AbstractConstructiveGeometry{T}, primitive, n_samples, CSG_scale, slice_val, plotattributes) where {T}
    projections = [:x, :y, :z]
    if plotattributes[:seriestype] in vcat(projections,:samplesurface)
        spacing = T((ismissing(CSG_scale) ? get_scale(csg) : CSG_scale)/n_samples)
        samples_csg = filter(p -> in(p,csg), sample(primitive, spacing))
        if plotattributes[:seriestype] in projections
            samples = filter(pt -> abs(getproperty(pt, plotattributes[:seriestype]) - slice_val) < spacing/2, samples_csg)
            proj = filter(x -> x != plotattributes[:seriestype], projections)
            internal_length_unit*getproperty.(samples, proj[1]), internal_length_unit*getproperty.(samples, proj[2])
        else
            samples_csg
        end
    else
        primitive
    end
end

@recipe function f(csg::AbstractConstructiveGeometry{T}; n_samples = 40, CSG_scale = missing, slice_val = T(0)) where {T}
    seriestype --> :csg
    unitformat --> :slash
    projections = [:x, :y, :z]
    ps = primitives(csg)
    @series begin
        label --> "CSG"
        if plotattributes[:seriestype] in projections
            proj = filter(x -> x != plotattributes[:seriestype], projections)
            xunit --> internal_length_unit
            yunit --> internal_length_unit
            xguide --> string(proj[1])
            yguide --> string(proj[2])
        else
            if !isClosedPrimitive(ps[1])
                fillcolor := :white
                fillalpha --> 0.2
                linewidth := 1
            end
        end
        _get_csg_plot_objects(csg, ps[1], n_samples, CSG_scale, slice_val, plotattributes)
    end
    for i in 2:length(ps)
        @series begin
            label := ""
            if plotattributes[:seriestype] in projections
                proj = filter(x -> x != plotattributes[:seriestype], projections)
                xunit --> internal_length_unit
                yunit --> internal_length_unit
                xguide --> string(proj[1])
                yguide --> string(proj[2])
            else
                if !isClosedPrimitive(ps[i])
                    fillcolor := :white
                    fillalpha --> 0.2
                    linewidth := 1
                end
            end
            _get_csg_plot_objects(csg, ps[i], n_samples, CSG_scale, slice_val, plotattributes)
        end
    end
end