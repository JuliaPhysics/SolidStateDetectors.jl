# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors
using Unitful
using StaticArrays
using LinearAlgebra

T = Float32

@testset "Shell Structure" begin
    center = CartesianPoint{T}(0, 0, 0)

    structure = ((4, T(0), T(0)),)
    pts = SolidStateDetectors.get_shell_structure_locations(structure, center, T(1))

    @test pts isa SVector{4,CartesianPoint{T}}
    @test length(pts) == 4

    # All points lie on sphere of radius 1
    @test all(abs(norm(p - center) - 1) < 1e-5 for p in pts)

    # Center of mass is at the center
    v = sum((p - center) for p in pts)
    com = center + v / T(4)
    @test all(abs(com[i] - center[i]) < 1e-5 for i in 1:3)
end

function test_platonic(constructor, N)
    center = CartesianPoint{T}(0,0,0)
    length_val = T(1.7)

    solid = constructor(center, length_val)

    @test solid isa SolidStateDetectors.PlatonicSolid{N,T}
    @test length(solid.locations) == N

    # Radius check
    @test all(abs(norm(p - center) - length_val) < 1e-4 for p in solid.locations)

    # Center of mass
    v = sum((p - center) for p in solid.locations)
    com = center + v / T(N)
    @test all(abs(com[i] - center[i]) < 1e-4 for i in 1:3)

    return true
end

@testset "Tetrahedron" begin
    test_platonic(SolidStateDetectors.Tetrahedron, 4)
end

@testset "Hexahedron" begin
    test_platonic(SolidStateDetectors.Hexahedron, 8)
end

@testset "Octahedron" begin
    test_platonic(SolidStateDetectors.Octahedron, 6)
end

@testset "Icosahedron" begin
    test_platonic(SolidStateDetectors.Icosahedron, 12)
end

@testset "Dodecahedron" begin
    test_platonic(SolidStateDetectors.Dodecahedron, 20)
end

@testset "Vertices" begin
    @test SolidStateDetectors.get_vertices(SolidStateDetectors.PlatonicSolid{4,T}) == 4
    @test SolidStateDetectors.get_vertices(SolidStateDetectors.PlatonicSolid{20,T}) == 20
end

@testset "PointCharge" begin
    c = CartesianPoint{T}(1,2,3)

    pc1 = SolidStateDetectors.PointCharge(c)
    @test pc1 isa SolidStateDetectors.PointCharge{T}
    @test length(pc1.locations) == 1
    @test pc1.locations[1] == c

    pc2 = SolidStateDetectors.PointCharge([c])
    @test pc2 isa SolidStateDetectors.PointCharge{T}
    @test length(pc2.locations) == 1
    @test pc2.locations[1] == c

    @test norm(pc1.locations[1] - c) == 0
end

@testset "Distance Line" begin
    # Line along the x-axis, origin at (0,0,0)
    a = CartesianPoint{T}(0,0,0)
    b = CartesianPoint{T}(1,0,0)
    l = SolidStateDetectors.ConstructiveSolidGeometry.Line(a, b)

    # Point on the line
    pt_on_line = CartesianPoint{T}(2,0,0)
    @test SolidStateDetectors.ConstructiveSolidGeometry.distance(pt_on_line, l) ≈ T(0)

    # Point above the line in z direction
    pt_above = CartesianPoint{T}(2,0,1)
    @test SolidStateDetectors.ConstructiveSolidGeometry.distance(pt_above, l) ≈ T(1)

    # Point above the line in y direction
    pt_side = CartesianPoint{T}(3,2,0)
    @test SolidStateDetectors.ConstructiveSolidGeometry.distance(pt_side, l) ≈ T(2)

    # Line along arbitrary direction
    a2 = CartesianPoint{T}(1,1,1)
    b2 = CartesianPoint{T}(2,3,2)
    l2 = SolidStateDetectors.ConstructiveSolidGeometry.Line(a2, b2)
    pt2 = CartesianPoint{T}(3,1,1)
    # Compute expected distance manually
    v = pt2 - l2.origin
    d_expected = norm(cross(v, l2.direction)) / norm(l2.direction)
    @test SolidStateDetectors.ConstructiveSolidGeometry.distance(pt2, l2) ≈ d_expected
end

@testset "Distance Edge" begin
    # Edge along x-axis
    a = CartesianPoint{T}(0,0,0)
    b = CartesianPoint{T}(1,0,0)
    e = SolidStateDetectors.ConstructiveSolidGeometry.Edge(a, b)

    # Point on the segment
    pt_on = CartesianPoint{T}(0.5,0,0)
    @test SolidStateDetectors.ConstructiveSolidGeometry.distance(pt_on, e) ≈ T(0)

    # Point beyond b: distance = ||pt - b||
    pt_beyond = CartesianPoint{T}(2,1,0)
    @test SolidStateDetectors.ConstructiveSolidGeometry.distance(pt_beyond, e) ≈ norm(pt_beyond - b)

    # Point before a: distance = ||pt - a||
    pt_before = CartesianPoint{T}(-1,1,0)
    @test SolidStateDetectors.ConstructiveSolidGeometry.distance(pt_before, e) ≈ norm(pt_before - a)

    # Point off the segment but projection falls on segment
    pt_off = CartesianPoint{T}(0.5,2,0)
    # perpendicular distance to line
    v = b - a
    expected = norm(cross(pt_off - a, v)) / norm(v)
    @test SolidStateDetectors.ConstructiveSolidGeometry.distance(pt_off, e) ≈ expected

    # Z axis edge
    e_vert = SolidStateDetectors.ConstructiveSolidGeometry.Edge(CartesianPoint{T}(0,0,0), CartesianPoint{T}(0,0,1))
    pt_off_vert = CartesianPoint{T}(1,0,0.5)
    @test SolidStateDetectors.ConstructiveSolidGeometry.distance(pt_off_vert, e_vert) ≈ T(1)
end
