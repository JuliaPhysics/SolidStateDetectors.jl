# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors: CartesianPoint, CartesianVector, CylindricalPoint, CylindricalVector
using StaticArrays: Size, SVector, SMatrix

@testset "points_and_vectors" begin
    a = CartesianPoint(1.0, 2.0, 3.0)
    b = CartesianPoint(3.0, 1.0, 2.0)
    v = CartesianVector(0.1, 0.2, 0.3)
    A = SMatrix{3,3}(1,2,3,4,5,6,7,8,9)

    @test @inferred(a + v) == CartesianPoint(1.1, 2.2, 3.3)
    @test @inferred(a - v) == CartesianPoint(0.9, 1.8, 2.7)
    @test @inferred(a - v) == CartesianPoint(0.9, 1.8, 2.7)
    @test @inferred(a - b) == CartesianVector(-2.0, 1.0, 1.0)
    @test @inferred(A * a) == CartesianPoint(30.0, 36.0, 42.0)

    @test @inferred(Size(a)) === Size(v)
    @test @inferred(size(a)) === size(v)
    @test @inferred(length(a)) === length(v)
    @test @inferred(CartesianPoint(a[1], a[2], a[3])) === a
    @test @inferred(CartesianPoint(a[1], a[2], a[3])) == a
    @test @inferred(CartesianPoint(a[1], a[2], a[3])) ≈ a


    a = CylindricalPoint(1.0, 2.0, 3.0)
    b = CylindricalPoint(3.0, 1.0, 2.0)

    v = CylindricalVector(0.1, 0.2, 0.3)

    @test  @inferred(a + v) == CylindricalPoint(1.1, 2.2, 3.3)
    @test  @inferred(a - v) == CylindricalPoint(0.9, 1.8, 2.7)
    @test  @inferred(a - b) == CylindricalVector(-2.0, 1.0, 1.0)
    @test  @inferred(Size(a)) === Size(v)
    @test  @inferred(size(a)) === size(v)
    @test  @inferred(length(a)) === length(v)
    @test  @inferred(CylindricalPoint(a[1], a[2], a[3])) === a
    @test  @inferred(CylindricalPoint(a[1], a[2], a[3])) == a
    @test  @inferred(CylindricalPoint(a[1], a[2], a[3])) ≈ a
end
