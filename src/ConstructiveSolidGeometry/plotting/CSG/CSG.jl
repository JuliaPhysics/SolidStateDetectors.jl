@recipe function f(csg::ConstructiveSolidGeometry.CSGUnion)
    linecolor --> :black
    csguniontype --> :union
    @series begin 
        linestyle --> :solid
        csg.a
    end
    @series begin 
        linestyle --> :solid
        csg.b
    end
end
@recipe function f(csg::ConstructiveSolidGeometry.CSGDifference)
    linecolor --> :black
    @series begin
        linestyle --> :solid
        csguniontype --> :difference
        csg.a
    end
    @series begin
        linestyle --> :dash
        csg.b
    end
end
@recipe function f(csg::ConstructiveSolidGeometry.CSGIntersection)
    linecolor --> :black
    csguniontype --> :intersection
    @series begin
        linestyle --> :dot
        csguniontype --> :difference
        csg.a
    end
    @series begin
        linestyle --> :dot
        csg.b
    end
end