primitives(vp::AbstractVolumePrimitive) = (vp,)
function primitives(csg::AbstractConstructiveGeometry)
    ps = []
    push!(ps, primitives(csg.a)...)
    push!(ps, primitives(csg.b)...)
    ps
end

@recipe function f(csg::AbstractConstructiveGeometry{T}; n_samples = 40) where {T}
    seriestype --> :csg
    ps = primitives(csg)
    spacing = (haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :samplesurface) ?  T(get_scale(csg)/n_samples) : nothing
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
                            linewidth := 0.5
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
                                linewidth := 0.5
                            end
                        end
                    end
                    ps[i]
                end
            end
        end
    end
end

