# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors.ConstructiveSolidGeometry: CartesianPoint, CartesianVector, 
    CartesianZero, cartesian_zero, CylindricalPoint, LocalAffineFrame, global_frame, frame_transformation, barycenter, geom_round
using StaticArrays: Size, SVector, SMatrix
using InverseFunctions: inverse

using Unitful

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

        @test @inferred(zero(a) + (a - zero(a))) == a

        @test @inferred(Size(a)) === Size(v)
        @test @inferred(size(a)) === size(v)
        @test @inferred(length(a)) === length(v)
        @test @inferred(CartesianPoint(a[1], a[2], a[3])) === a
        @test @inferred(CartesianPoint(a[1], a[2], a[3])) == a
        @test @inferred(CartesianPoint(a[1], a[2], a[3])) ≈ a

        frame = LocalAffineFrame(b, A)

        f = frame_transformation(frame, global_frame)
        @test @inferred(f(a)) == cartesian_zero + A * (a - cartesian_zero) + CartesianVector(b[1], b[2], b[3])
        @test @inferred(inverse(f)(f(a))) ≈ a


        #@test @inferred(CartesianZero{Float32}() * u"mm") === CartesianZero{typeof(zero(Float32) * u"mm")}()

        A = [CartesianPoint{Float32}(x,0,0) for x in -2:2]
        @test isapprox(barycenter(A), CartesianPoint{Float32}(0,0,0))

        S = SVector{length(A)}(A)
        @test isapprox(barycenter(S), CartesianPoint{Float32}(0,0,0))

        # test types and units
        @test CartesianPoint(1, 2, 3) isa CartesianPoint{Float64}
        @test CartesianPoint(1, 2, 3f0) isa CartesianPoint{Float32}
        @test CartesianPoint(1, 2.0f0, Float16(3)) isa CartesianPoint{Float32}
        @test CartesianPoint(1.0, 2.0f0, Float16(3)) isa CartesianPoint{Float64}
        @test CartesianPoint(1, 2, Float16(3)) isa CartesianPoint{Float16}
        @test CartesianPoint(1.0u"m", 2.0u"m", 3.0f0u"m") isa CartesianPoint{Float64}
        @test CartesianPoint(1.0u"m", 2.0f0u"m", 3.0f0u"m") isa CartesianPoint{Float64}
        @test CartesianPoint(1.0f0u"m", 2.0f0u"m", 3.0f0u"m") isa CartesianPoint{Float32}
        @test CartesianPoint(1u"m", 2u"m", 3u"m") isa CartesianPoint{Float64}
        @test CartesianPoint(1.0u"mm", 2.0u"cm", 3.0f0u"m") isa CartesianPoint{Float64}
        @test CartesianPoint(1.0f0u"mm", 2.0f0u"mm", 3.0u"m") isa CartesianPoint{Float64}
        @test CartesianPoint(1.0u"mm", 2.0f0u"cm", 3.0f0u"m") isa CartesianPoint{Float64}
        @test CartesianPoint(1.0f0u"mm", 2.0f0u"cm", 3.0f0u"m") isa CartesianPoint{Float32}
        @test CartesianPoint(1u"mm", 2u"cm", 3u"m") isa CartesianPoint{Float64}
        @test CartesianPoint(1u"mm", 0, 0) isa CartesianPoint{Float64}
        @test CartesianPoint(1u"m", 2u"m", 3f0u"m") isa CartesianPoint{Float32}
        
        # test throwing errors with wrong units
        @test_throws ArgumentError CartesianPoint(1u"m", 2u"rad", 3u"m")
        @test_throws ArgumentError CartesianPoint(1u"s", 2u"m", 3u"m")
        @test_throws ArgumentError CartesianPoint(1u"m", 2u"m", 3u"kg")


        # test Base.convert
        p32 = CartesianPoint{Float32}(1.0, 2.0, 3.0)
        p64 = CartesianPoint{Float64}(1.0, 2.0, 3.0)
        @test @inferred(convert(CartesianPoint{Float64}, p32)) isa CartesianPoint{Float64}
        @test @inferred(convert(CartesianPoint{Float64}, p32)) == p64

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

        A = [CylindricalPoint{Float32}(x,0,0) for x in -2:2]
        @test isapprox(barycenter(A), CylindricalPoint{Float32}(0,0,0))

        S = SVector{length(A)}(A)
        @test isapprox(barycenter(S), CylindricalPoint{Float32}(0,0,0))

        # test types and units
        @test CylindricalPoint(1, 2, 3) isa CylindricalPoint{Float64}
        @test CylindricalPoint(1, 2, 3f0) isa CylindricalPoint{Float32}
        @test CylindricalPoint(1, 2.0f0, Float16(3)) isa CylindricalPoint{Float32}
        @test CylindricalPoint(1.0, 2.0f0, Float16(3)) isa CylindricalPoint{Float64}
        @test CylindricalPoint(1, 2, Float16(3)) isa CylindricalPoint{Float16}
        @test CylindricalPoint(1.0u"m", 2.0u"rad", 3.0f0u"m") isa CylindricalPoint{Float64}
        @test CylindricalPoint(1.0u"m", 2.0f0u"rad", 3.0f0u"m") isa CylindricalPoint{Float64}
        @test CylindricalPoint(1.0f0u"m", 2.0f0u"rad", 3.0f0u"m") isa CylindricalPoint{Float32}
        @test CylindricalPoint(1u"m", 2u"rad", 3u"m") isa CylindricalPoint{Float64}
        @test CylindricalPoint(1.0u"mm", 2.0u"rad", 3.0f0u"m") isa CylindricalPoint{Float64}
        @test CylindricalPoint(1.0f0u"mm", 2.0f0u"rad", 3.0u"m") isa CylindricalPoint{Float64}
        @test CylindricalPoint(1.0u"mm", 2.0f0u"rad", 3.0f0u"m") isa CylindricalPoint{Float64}
        @test CylindricalPoint(1.0f0u"mm", 2.0f0u"rad", 3.0f0u"m") isa CylindricalPoint{Float32}
        @test CylindricalPoint(1u"mm", 2u"rad", 3u"m") isa CylindricalPoint{Float64}
        @test CylindricalPoint(1u"mm", 0, 0) isa CylindricalPoint{Float64}
        @test CylindricalPoint(1u"m", 2u"rad", 3f0u"m") isa CylindricalPoint{Float32}
        
        # test throwing errors with wrong units
        @test_throws ArgumentError CylindricalPoint(1u"rad", 2u"rad", 3u"m")
        @test_throws ArgumentError CylindricalPoint(1u"m", 2u"m", 3u"m")
        @test_throws ArgumentError CylindricalPoint(1u"m", 2u"rad", 3u"rad")
    end
end

@testset "CartesianVector" begin

    @test CartesianVector(1, 2, 3) isa CartesianVector{Float64}
    @test CartesianVector(1, 2, 3f0) isa CartesianVector{Float32}
    @test CartesianVector(1, 2.0f0, Float16(3)) isa CartesianVector{Float32}
    @test CartesianVector(1.0, 2.0f0, Float16(3)) isa CartesianVector{Float64}
    @test CartesianVector(1, 2, Float16(3)) isa CartesianVector{Float16}

    @test CartesianVector(1.0u"m", 2.0u"m", 3.0f0u"m") isa CartesianVector{Float64}
    @test CartesianVector(1.0u"m", 2.0f0u"m", 3.0f0u"m") isa CartesianVector{Float64}
    @test CartesianVector(1.0f0u"m", 2.0f0u"m", 3.0f0u"m") isa CartesianVector{Float32}
    @test CartesianVector(1u"m", 2u"m", 3u"m") isa CartesianVector{Float64}

    @test CartesianVector(1.0u"mm", 2.0u"cm", 3.0f0u"m") isa CartesianVector{Float64}
    @test CartesianVector(1.0f0u"mm", 2.0f0u"mm", 3.0u"m") isa CartesianVector{Float64}
    @test CartesianVector(1.0u"mm", 2.0f0u"cm", 3.0f0u"m") isa CartesianVector{Float64}
    @test CartesianVector(1.0f0u"mm", 2.0f0u"cm", 3.0f0u"m") isa CartesianVector{Float32}
    @test CartesianVector(1u"mm", 2u"cm", 3u"m") isa CartesianVector{Float64}

    @test CartesianVector(1u"mm", 0, 0) isa CartesianVector{Float64}
    @test CartesianVector(1u"m", 2u"m", 3f0u"m") isa CartesianVector{Float32}

    @test_throws ArgumentError CartesianVector(1u"m", 2u"rad", 3u"m")
    @test_throws ArgumentError CartesianVector(1u"s", 2u"m", 3u"m")
    @test_throws ArgumentError CartesianVector(1u"m", 2u"m", 3u"kg")

    v1 = CartesianVector(1, 2, 3)
    @test v1 == CartesianVector(1.0, 2.0, 3.0)

    v2 = CartesianVector{Float32}(x=1, y=2, z=3)
    @test v2 == CartesianVector(1f0, 2f0, 3f0)

    z = zero(CartesianVector{Float64})
    @test z == CartesianVector(0.0, 0.0, 0.0)

    v = CartesianVector(1.0, 2.0, 3.0)
    v_mul = v * u"m"
    @test v_mul == CartesianVector(1.0u"m", 2.0u"m", 3.0u"m")

    v = CartesianVector(1f0, 2f0, 3f0)
    @test -v == CartesianVector(-1f0, -2f0, -3f0)

    cz1 = CartesianZero{Float32}()
    cz2 = CartesianZero{Float64}()
    @test -(cz1, cz2) == CartesianVector(0.0, 0.0, 0.0)
end

@testset "affine operations" begin
    # Test affine geometrical operations
    pt = CartesianPoint(1.0, 2.0, 3.0)
    v = CartesianVector(0.1, -0.2, 0.3)

    @test pt - cartesian_zero isa CartesianVector
    @test cartesian_zero + v isa CartesianPoint
    @test pt - v isa CartesianPoint
    @test pt + v isa CartesianPoint
    @test transpose(pt) isa CartesianPoint
    @test adjoint(pt) isa CartesianPoint
    @test CartesianPoint{Float32}(pt) isa CartesianPoint{Float32}
    @test Base.convert(CartesianPoint{Float32}, pt) isa CartesianPoint{Float32}
    @test_throws DimensionMismatch pt + [1.0, 2.0]
    @test_throws DimensionMismatch pt - [1.0, 2.0]

    # Test frame_transformation accepts AbstractVector{<:AbstractCoordinateVector}
    vectors = [CartesianVector(1.0, 0.0, 0.0), CartesianVector(0.0, 1.0, 0.0)]
    transformed = frame_transformation(LocalAffineFrame(cartesian_zero, one(SMatrix{3,3,Float64,9})), global_frame, vectors)
    @test transformed isa AbstractVector{<:CartesianVector}
    @test transformed == vectors 
end

@testset "geom_round" begin
    for T in (Float32, Float64)
        pts = [CartesianPoint{T}(1.0, 1e-13, 2.0), CylindricalPoint{T}(1.0, 1e-13, 0.0)]
        rounded = geom_round(pts)
        @test all(isa.(rounded, typeof.(pts)))
        @test all(rounded .!= pts)
        @test all(iszero.(getindex.(rounded, 2)))
    end
end
