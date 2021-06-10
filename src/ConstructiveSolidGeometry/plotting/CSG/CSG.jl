@recipe function f(csg::ConstructiveSolidGeometry.CSGUnion)
    linecolor --> :black
    csguniontype --> :union
    csg.a, csg.b
end
@recipe function f(csg::ConstructiveSolidGeometry.CSGDifference)
    linecolor --> :black
    csguniontype --> :difference
    csg.a, csg.b
end
@recipe function f(csg::ConstructiveSolidGeometry.CSGIntersection)
    linecolor --> :black
    csguniontype --> :intersection
    csg.a, csg.b
end