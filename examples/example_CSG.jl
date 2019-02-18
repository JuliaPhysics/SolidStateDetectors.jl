using Test
using BenchmarkTools, LinearAlgebra, StaticArrays
using Plots; pyplot()
using SolidStateDetectors

T = Float64

np = SSD.CartesianPoint{T}(0, 0, 0)
pt = SSD.CartesianPoint{T}(0.64, 0.15, 0.54)

ex = SSD.CartesianVector{T}(1.0, 0.0, 0.0)
ey = SSD.CartesianVector{T}(0.0, 1.0, 0.0)
ez = SSD.CartesianVector{T}(0.0, 0.0, 1.0)

rect = SSD.RectangularCuboid{T}( np, ex, ey, ez )
@test in(np, rect)
@test in(np + ex, rect)
@test in(np + ey, rect)
@test in(np + ez, rect)
@test in(pt, rect)
@test !in(pt + ex, rect)
@test !in(pt + ey, rect)
@test !in(pt + ez, rect)


cyl = SSD.Cylinder{T}( np, ez, T(0.4) )


l1 = SSD.Line( SSD.CartesianPoint{T}(1, 0, 0), SSD.CartesianVector{T}(0, 1, 0) )
l2 = SSD.Line( SSD.CartesianPoint{T}(0, 1, 0), SSD.CartesianVector{T}(1, 0, 0) )
r1 = SSD.Ray( np, ex )
r2 = SSD.Ray( np, ey )
edge = SSD.LineSegment(np, SSD.CartesianPoint{T}(0, 1, 0))
edge1 = SSD.LineSegment(SSD.CartesianPoint{T}(0, 1, 0),  SSD.CartesianPoint{T}(1, 0, 0))
edge2 = SSD.LineSegment(SSD.CartesianPoint{T}(-10, 0, 0),  SSD.CartesianPoint{T}(1, 1, 0))


csg_diff = rect - cyl
@test !in(np, csg_diff)
@test in(np + ey, csg_diff)