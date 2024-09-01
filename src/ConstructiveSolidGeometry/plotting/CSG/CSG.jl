@inline primitives(vp::AbstractVolumePrimitive) = (vp,)
@inline function primitives(csg::AbstractConstructiveGeometry)
    (
        primitives(csg.a)...,
        primitives(csg.b)...
    )
end

@recipe function f(csg::AbstractConstructiveGeometry{T}; n_samples = 40, CSG_scale = missing, slice_val = T(0)) where {T}
    seriestype --> :csg
    xunit --> internal_length_unit
    yunit --> internal_length_unit
    projections = [:x, :y, :z]
    ps = primitives(csg)
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] in vcat(projections,:samplesurface)
        spacing::T = T((ismissing(CSG_scale) ? get_scale(csg) : CSG_scale)/n_samples)
    end
    @series begin
        label --> "CSG"
        if haskey(plotattributes, :seriestype) 
            if plotattributes[:seriestype] == :samplesurface
                zunit --> internal_length_unit
                filter(p -> in(p,csg), sample(ps[1], spacing))
            elseif plotattributes[:seriestype] in projections
                #Unlike f(CartesianPoint), there is no ssd plot recipe for f(array,array) so must define attributes here
                xguide --> "x"
                yguide --> "y"
                unitformat --> :slash
                samples = filter(pt -> abs(getproperty(pt, plotattributes[:seriestype]) - slice_val) < spacing/2, filter(p -> in(p,csg), sample(ps[1], spacing)))
                proj = filter(x -> x != plotattributes[:seriestype], projections)
                internal_length_unit*getproperty.(samples, proj[1]), internal_length_unit*getproperty.(samples, proj[2])
            else
                zunit --> internal_length_unit
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
                elseif plotattributes[:seriestype] in projections
                    xguide --> "x"
                    yguide --> "y"
                    unitformat --> :slash
                    samples = filter(pt -> abs(getproperty(pt, plotattributes[:seriestype]) - slice_val) < spacing/2, filter(p -> in(p,csg), sample(ps[i], spacing)))
                    proj = filter(x -> x != plotattributes[:seriestype], projections)
                    internal_length_unit*getproperty.(samples, proj[1]), internal_length_unit*getproperty.(samples, proj[2])
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

