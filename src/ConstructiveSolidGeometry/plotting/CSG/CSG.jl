primitives(vp::AbstractVolumePrimitive) = (vp,)
function primitives(csg::AbstractConstructiveGeometry)
    ps = []
    push!(ps, primitives(csg.a)...)
    push!(ps, primitives(csg.b)...)
    ps
end

@recipe function f(csg::AbstractConstructiveGeometry)
    ps = primitives(csg)
    @series begin
        label --> "CSG"
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
    for i in 2:length(ps)
        @series begin
            label := ""
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

