primitives(vp::AbstractVolumePrimitive) = (vp,)
function primitives(csg::AbstractConstructiveGeometry)
    ps = []
    push!(ps, primitives(csg.a)...)
    push!(ps, primitives(csg.b)...)
    ps
end

extremum(csg::AbstractConstructiveGeometry{T}) where {T} = maximum([extremum(s) + norm(s.origin) for p in primitives(csg) for s in surfaces(p)])

@recipe function f(csg::AbstractConstructiveGeometry{T}; n_samples = 200) where {T}
    ps = primitives(csg)
    spacing = (haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :samplesurface) ?  T(extremum(csg)/n_samples) : nothing
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
        if haskey(plotattributes, :seriestype) 
            if plotattributes[:seriestype] == :samplesurface
                seriestype := :scatter
                seriescolor --> 1
                seriesalpha --> 0.2
                markerstrokewidth --> 0
                filter(p -> in(p,csg,csgtol = 10000*csg_default_tol(T)), vertices(ps[1], spacing))
            else
                ps[1]
            end
        end
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
            if haskey(plotattributes, :seriestype) 
                if plotattributes[:seriestype] == :samplesurface
                    seriestype := :scatter
                    seriescolor --> 1
                    seriesalpha --> 0.2
                    markerstrokewidth --> 0
                    filter(p -> in(p,csg,csgtol = 10000*csg_default_tol(T)), vertices(ps[i], spacing))
                else
                    ps[i]
                end
            end
        end
    end
end

