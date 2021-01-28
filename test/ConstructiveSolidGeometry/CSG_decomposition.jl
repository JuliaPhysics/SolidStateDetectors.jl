using Test
using SolidStateDetectors
using LinearAlgebra
using StaticArrays
using Rotations

import SolidStateDetectors.ConstructiveSolidGeometry: Geometry
using SolidStateDetectors.ConstructiveSolidGeometry: 
    Box, Cone, HexagonalPrism, Sphere, Tube,
    CartesianVector, scale, get_decomposed_volumes

T = Float64

example_primitive_dir = joinpath(@__DIR__, "../../examples/example_primitive_files")
@testset "Test decomposition" begin
    rot_x = RotMatrix(RotX(deg2rad(90)))
    rot_x_inv = inv(rot_x)
    scaling = SVector{3,T}(1,2,0.5)
    t = CartesianVector{T}(1,0,0)
    t_inv = -t
    tube = Tube(T(1),T(2)) # r from 0..1, z from -1..1
    sphere = Sphere(T(1)) # r from 0..1
    box = Box(T(1), T(2), T(3)) # x from -0.5..0.5, y from -1.0..1.0, z from -1.5..1.5
    @testset "Combinations with only positive geometries" begin
        g1 = box + (sphere & tube)
        pos, neg = get_decomposed_volumes(g1)
        @test box in pos
        @test sphere in pos
        @test tube in pos
        @test !(box in neg)
        @test !(sphere in neg)
        @test !(tube in neg)
        g2 = (box + t) + ((sphere + t_inv) & (rot_x * tube))
        pos, neg = get_decomposed_volumes(g2)
        @test (box + t) in pos
        @test (sphere + t_inv) in pos
        @test (rot_x * tube) in pos
        @test !((box + t) in neg)
        @test !((sphere + t_inv) in neg)
        @test !((rot_x * tube) in neg)
    end
    @testset "Combinations with positive and negative geometries" begin
        g1 = box + sphere - tube
        pos, neg = get_decomposed_volumes(g1)
        @test box in pos
        @test sphere in pos
        @test !(tube in pos)
        @test !(box in neg)
        @test !(sphere in neg)
        @test tube in neg
        g2 = (box + t) - ((sphere + t_inv) & (rot_x * tube))
        pos, neg = get_decomposed_volumes(g2)
        @test (box + t) in pos
        @test !((sphere + t_inv) in pos)
        @test !((rot_x * tube) in pos)
        @test !((box + t) in neg)
        @test (sphere + t_inv) in neg
        @test (rot_x * tube) in neg
    end
    @testset "Combined rotation of set" begin
        g1 = rot_x * (box + sphere - tube)
        pos, neg = get_decomposed_volumes(g1)
        @test rot_x * box in pos
        @test rot_x * sphere in pos
        @test !(rot_x * tube in pos)
        @test !(rot_x * box in neg)
        @test !(rot_x * sphere in neg)
        @test rot_x * tube in neg
        g2 = scale(rot_x * (((box + t) - ((sphere + t_inv) & (rot_x * tube))) + t), scaling)
        pos, neg = get_decomposed_volumes(g2)
        @test scale(rot_x * (box + 2t), scaling) in pos
        @test !(scale(rot_x * sphere, scaling) in pos)
        @test !(scale(rot_x * (rot_x * tube + t), scaling) in pos)
        @test !(scale(rot_x * (box + 2t), scaling) in neg)
        @test scale(rot_x * sphere, scaling) in neg
        @test scale(rot_x * (rot_x * tube + t), scaling) in neg
    end
end