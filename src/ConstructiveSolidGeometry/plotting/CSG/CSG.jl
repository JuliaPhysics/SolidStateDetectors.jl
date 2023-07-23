primitives(vp::AbstractVolumePrimitive) = (vp,)
function primitives(csg::AbstractConstructiveGeometry)
    (
        (@inline primitives(csg.a))...,
        (@inline primitives(csg.b))...
    )
end

@recipe function f(csg::AbstractConstructiveGeometry{T}; n_samples = 40, CSG_scale = missing) where {T}
    seriestype --> :csg
    xunit --> internal_length_unit
    yunit --> internal_length_unit
    zunit --> internal_length_unit
    ps = primitives(csg)
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :samplesurface
        spacing::T = T((ismissing(CSG_scale) ? get_scale(csg) : CSG_scale)/n_samples)
    end
    @series begin
        label --> "CSG"
        if haskey(plotattributes, :seriestype) 
            if plotattributes[:seriestype] == :samplesurface
                if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
                    aspect_ratio --> 1.0
                end 
                seriesalpha --> 0.2
                filter(p -> in(p,csg), sample(ps[1], spacing))
            else
                if !isClosedPrimitive(ps[1])
                    fillcolor := :white
                    fillalpha --> 0.2
                    if haskey(plotattributes, :seriestype) 
                        if plotattributes[:seriestype] in [:csg, :wireframe] 
                            linewidth := 1
                        end
                    end
                end
                ps[1]
            end
        end
    end
    for i in 2:length(ps)
        @series begin
            label := ""
            if haskey(plotattributes, :seriestype) 
                if plotattributes[:seriestype] == :samplesurface
                    seriesalpha --> 0.2
                    filter(p -> in(p,csg), sample(ps[i], spacing))
                else
                    if !isClosedPrimitive(ps[i])
                        fillcolor := :white
                        fillalpha --> 0.2
                        if haskey(plotattributes, :seriestype) 
                            if plotattributes[:seriestype] in [:csg, :wireframe] 
                                linewidth := 1
                            end
                        end
                    end
                    ps[i]
                end
            end
        end
    end
end

