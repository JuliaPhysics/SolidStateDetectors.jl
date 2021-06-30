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
        linestyle := isClosedPrimitive(ps[1]) ? :solid : :dot
        ps[1]
    end
    for i in 2:length(ps)
        @series begin
            label := nothing
            linestyle := isClosedPrimitive(ps[i]) ? :solid : :dot
            ps[i]
        end
    end
end

