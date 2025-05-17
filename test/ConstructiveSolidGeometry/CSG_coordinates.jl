# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors: ConstructiveSolidGeometry as CSG
using .CSG: CartesianPoint, CartesianVector, CylindricalPoint
using .CSG: LocalAffineFrame, global_frame, frame_transformation
using StaticArrays: Size, SVector, SMatrix
using InverseFunctions: inverse

@testset "points_and_vectors" begin
    @testset "cartesian" begin
        cart = @inferred CartesianVector(x=2f0,z=1f0)
        @inferred CartesianVector{Float32}(x=2)
        @test cart.x == Float32(2)
        cart = @inferred CartesianPoint(x=2f0,z=1f0)
        @inferred CartesianPoint{Float32}(x=2)

        a = CartesianPoint(1.0, 2.0, 3.0)
        b = CartesianPoint(3.0, 1.0, 2.0)
        v = CartesianVector(0.1, 0.2, 0.3)
        A = SMatrix{3,3}(1.0, 4.0, 7.0, 4.0, 5.0, 8.0, 7.0, 8.0, 9.0)

        @test @inferred(a + v) == CartesianPoint(1.1, 2.2, 3.3)
        @test @inferred(a - v) == CartesianPoint(0.9, 1.8, 2.7)
        @test @inferred(a - v) == CartesianPoint(0.9, 1.8, 2.7)
        @test @inferred(a - b) == CartesianVector(-2.0, 1.0, 1.0)
        @test @inferred(A * a) == CartesianPoint(30.0, 38.0, 50.0)

        @test @inferred(zero(a) + (a - zero(a))) == a

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
    end

    @testset "cylindrical" begin
        cyl = @inferred CylindricalPoint{Float32}(r=2.,z=1.)
        cyl2 = @inferred CylindricalPoint(φ=3π)
        @test CartesianPoint(cyl) == CartesianPoint(x=2f0,z=1f0)

        a = CylindricalPoint(1.0, π/2, 3.0)
        b = CylindricalPoint(3.0, π/2, 2.0)
        v = CartesianVector(0.0, 0.2, 0.3)

        @test  @inferred(a + v) ≈ CylindricalPoint(1.2, π/2, 3.3)
        @test  @inferred(a - v) ≈ CylindricalPoint(0.8, π/2, 2.7)
        @test  @inferred(a - b) ≈ CartesianVector(0, -2, 1)

        @test @inferred(zero(a) + (a - zero(a))) == a

        @test  @inferred(Size(a)) === Size(v)
        @test  @inferred(size(a)) === size(v)
        @test  @inferred(length(a)) === length(v)
        @test  @inferred(CylindricalPoint(a[1], a[2], a[3])) === a
        @test  @inferred(CylindricalPoint(a[1], a[2], a[3])) == a
        @test  @inferred(CylindricalPoint(a[1], a[2], a[3])) ≈ a
    end
end
