# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors: CartesianPoint, CartesianVector, CylindricalPoint
using SolidStateDetectors.ConstructiveSolidGeometry: LocalAffineFrame, global_frame, frame_transformation
using StaticArrays: Size, SVector, SMatrix
using InverseFunctions: inverse

@testset "points_and_vectors" begin
    a = CartesianPoint(1.0, 2.0, 3.0)
    b = CartesianPoint(3.0, 1.0, 2.0)
    v = CartesianVector(0.1, 0.2, 0.3)
    A = SMatrix{3,3}(1.0, 4.0, 7.0, 4.0, 5.0, 8.0, 7.0, 8.0, 9.0)

    @test @inferred(a + v) == CartesianPoint(1.1, 2.2, 3.3)
    @test @inferred(a - v) == CartesianPoint(0.9, 1.8, 2.7)
    @test @inferred(a - v) == CartesianPoint(0.9, 1.8, 2.7)
    @test @inferred(a - b) == CartesianVector(-2.0, 1.0, 1.0)
    @test @inferred(A * a) == CartesianPoint(30.0, 38.0, 50.0)

    @test @inferred(Size(a)) === Size(v)
    @test @inferred(size(a)) === size(v)
    @test @inferred(length(a)) === length(v)
    @test @inferred(CartesianPoint(a[1], a[2], a[3])) === a
    @test @inferred(CartesianPoint(a[1], a[2], a[3])) == a
    @test @inferred(CartesianPoint(a[1], a[2], a[3])) ≈ a


    frame = LocalAffineFrame(b, A)

    f = frame_transformation(frame, global_frame)
    @test @inferred(f(a)) == A * a + CartesianVector(b[1], b[2], b[3])
    @test @inferred(inverse(f)(f(a))) ≈ a
    # ToDo: Add more tests!


    a = CylindricalPoint(1.0, 2.0, 3.0)
    b = CylindricalPoint(3.0, 1.0, 2.0)

    #@test  @inferred(a + v) == 
    #@test  @inferred(a - v) == 
    #@test  @inferred(a - b) == 
    @test  @inferred(Size(a)) === Size(v)
    @test  @inferred(size(a)) === size(v)
    @test  @inferred(length(a)) === length(v)
    @test  @inferred(CylindricalPoint(a[1], a[2], a[3])) === a
    @test  @inferred(CylindricalPoint(a[1], a[2], a[3])) == a
    @test  @inferred(CylindricalPoint(a[1], a[2], a[3])) ≈ a
end
