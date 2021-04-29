using Test
using SolidStateDetectors
using LinearAlgebra
using Rotations
using IntervalSets
using StaticArrays

using SolidStateDetectors.ConstructiveSolidGeometry:
    CartesianPoint, CartesianVector, CylindricalPoint, scale,
    _in_angular_interval_closed, _in_angular_interval_open,
    ConalPlane, ConeMantle, CylindricalAnnulus, ToroidalAnnulus, TorusMantle,
    Tube, Cone, Torus, Box, Sphere, HexagonalPrism

@testset "Test CSG" begin

    for T in [Float64, Float32]

        @testset  "Precision type: $(T)" begin

            @testset "Points" begin
                cp1 = CartesianPoint{T}(1, 2, 3)
                @test cp1 + cp1 == CartesianPoint{T}(2, 4, 6)
                @test cp1 - cp1 == CartesianPoint{T}(0, 0, 0)
                @test cp1 ⋅ cp1 == 14
                @test norm(cp1) == sqrt(T(14))
            end

            @testset "Angular Intervals" begin
                @testset "FullCircle" begin
                    full_2π = T(0)..T(2π)
                    @test _in_angular_interval_closed(T(0), full_2π)
                    @test _in_angular_interval_closed(T(2π), full_2π)
                    @test _in_angular_interval_closed(T(4π), full_2π)
                    @test _in_angular_interval_closed(T(-π), full_2π)
                    @test _in_angular_interval_closed(T(-2π), full_2π)
                    @test !_in_angular_interval_open(T(0), full_2π)
                    @test !_in_angular_interval_open(T(2π), full_2π)
                    @test !_in_angular_interval_open(T(4π), full_2π)
                    @test _in_angular_interval_open(T(-π), full_2π)
                    @test !_in_angular_interval_open(T(-2π), full_2π)
                end
                @testset "Partial Circle" begin
                    partial_2π = T(0)..T(π)
                    @test _in_angular_interval_closed(T(0), partial_2π)
                    @test _in_angular_interval_closed(T(2π), partial_2π)
                    @test _in_angular_interval_closed(T(π), partial_2π)
                    @test _in_angular_interval_closed(T(3π), partial_2π)
                    @test _in_angular_interval_closed(T(-π), partial_2π)
                    @test !_in_angular_interval_closed(T(-π/2), partial_2π)
                    @test _in_angular_interval_closed(T(π/2), partial_2π)
                    @test _in_angular_interval_closed(T(-3π/2), partial_2π)
                    @test !_in_angular_interval_open(T(0), partial_2π)
                    @test !_in_angular_interval_open(T(π), partial_2π)
                    #@test !_in_angular_interval_open(T(3π), partial_2π)
                    @test !_in_angular_interval_open(T(-π), partial_2π)
                    @test _in_angular_interval_open(T(π/2), partial_2π)
                    @test _in_angular_interval_open(T(-3π/2), partial_2π)
                end
                @testset "Negative Interval" begin
                    neg_interval = -T(π/2)..T(π/2)
                    @test _in_angular_interval_closed(T(0), neg_interval)
                    @test _in_angular_interval_closed(T(2π), neg_interval)
                    @test !_in_angular_interval_closed(T(π), neg_interval)
                    @test !_in_angular_interval_closed(T(3π), neg_interval)
                    @test !_in_angular_interval_closed(T(-π), neg_interval)
                    @test _in_angular_interval_closed(T(-π/2), neg_interval)
                    #@test _in_angular_interval_closed(T(-3π/2), neg_interval)
                    @test _in_angular_interval_closed(T(3π/2), neg_interval)
                    @test _in_angular_interval_closed(T(π/2), neg_interval)
                    @test _in_angular_interval_open(T(0), neg_interval)
                    @test !_in_angular_interval_open(T(π/2), neg_interval)
                    @test !_in_angular_interval_open(T(-π/2), neg_interval)
                    @test !_in_angular_interval_open(T(3π/2), neg_interval)
                    @test !_in_angular_interval_open(T(-3π/2), neg_interval)
                end
            end

            @testset "Surface Primitives" begin
                @testset "ConalPlane" begin
                    cplane = ConalPlane()
                    @test CartesianPoint{T}(0, 0, 0) in cplane
                    @test CartesianPoint{T}(1, 0, 0) in cplane
                    @test CartesianPoint{T}(1, 0, 0.5) in cplane
                    @test !(CartesianPoint{T}(1, 0.1, 0.5) in cplane)
                end
                @testset "ConeMantle" begin
                    cmantle = ConeMantle(T(0.0), T(2.0), T(2.0))
                    @test CartesianPoint{T}(0, 0, -1) in cmantle
                    @test CartesianPoint{T}(2, 0, 1) in cmantle
                    @test CartesianPoint{T}(0, 1, 0) in cmantle
                    @test !(CartesianPoint{T}(0, 0, 0) in cmantle)
                end
                @testset "CylindricalAnnulus" begin
                    cannulus = CylindricalAnnulus(T(2.0), T(0.0))
                    @test CartesianPoint{T}(0, 0, 0) in cannulus
                    @test CartesianPoint{T}(2, 0, 0) in cannulus
                    @test CartesianPoint{T}(0, -2, 0) in cannulus
                    @test !(CartesianPoint{T}(2, 2, 0) in cannulus)
                    @test !(CartesianPoint{T}(0, 0, 0.1) in cannulus)
                end
                @testset "ToroidalAnnulus" begin
                    tannulus = ToroidalAnnulus()
                    @test CartesianPoint{T}(0, 0, 0) in tannulus
                    @test CartesianPoint{T}(2, 0, 0) in tannulus
                    @test CartesianPoint{T}(1, 0, 1) in tannulus
                    @test CartesianPoint{T}(1, 0, -1) in tannulus
                    @test !(CartesianPoint{T}(0, 0.1, 0) in tannulus)
                    @test !(CartesianPoint{T}(2, 0, 0.1) in tannulus)
                end
                @testset "TorusMantle" begin
                    tmantle = TorusMantle()
                    @test CartesianPoint{T}(0, 0, 0) in tmantle
                    @test CartesianPoint{T}(2, 0, 0) in tmantle
                    @test CartesianPoint{T}(1, 0, 1) in tmantle
                    @test CartesianPoint{T}(1, 0, -1) in tmantle
                    @test CartesianPoint{T}(0, 2, 0) in tmantle
                    @test CartesianPoint{T}(0, 1, 1) in tmantle
                    @test CartesianPoint{T}(0, 1, -1) in tmantle
                    @test !(CartesianPoint{T}(0, 0.1, 0) in tmantle)
                    @test !(CartesianPoint{T}(2, 0, 0.1) in tmantle)
                end
            end

            @testset "Volume Primitives" begin
                @testset "Tube" begin
                    tube = Tube(T(0.2), T(2.0)) # r from 0..0.2, z from -1.0..1.0
                    @test CartesianPoint{T}(0, 0, 0) in tube
                    @test CartesianPoint{T}(tube.r, 0, tube.z) in tube
                    @test !(CartesianPoint{T}(tube.r, 0, 1.1*tube.z) in tube)
                    @test !(CartesianPoint{T}(1.1*tube.r, 0, tube.z) in tube)
                end
                @testset "Torus" begin
                    torus = Torus(T(2.0), T(1.0), T(0.0)) # r_torus = 2, r_tube = 1, z = 0
                    @test CartesianPoint{T}(torus.r_torus, 0, 0) in torus
                    @test CartesianPoint{T}(torus.r_torus + torus.r_tube, 0, 0) in torus
                    @test CartesianPoint{T}(torus.r_torus - torus.r_tube, 0, 0) in torus
                    @test CartesianPoint{T}(0, torus.r_torus + torus.r_tube, 0) in torus
                    @test CartesianPoint{T}(0, torus.r_torus - torus.r_tube, 0) in torus
                    @test CartesianPoint{T}(torus.r_torus, 0, torus.r_tube) in torus
                    @test CartesianPoint{T}(torus.r_torus, 0, -torus.r_tube) in torus
                    @test CartesianPoint{T}(0, torus.r_torus, torus.r_tube) in torus
                    @test CartesianPoint{T}(0, torus.r_torus, -torus.r_tube) in torus
                    @test !(CartesianPoint{T}(torus.r_torus + torus.r_tube, 0, 0.1) in torus)
                    @test !(CartesianPoint{T}(torus.r_torus + torus.r_tube + 0.1, 0, 0) in torus)
                end
                @testset "Box" begin
                    box = Box(T(1.0), T(2.0), T(3.0)) # x from -0.5..0.5, y from -1.0..1.0, z from -1.5..1.5
                    @test CartesianPoint{T}(0, 0, 0) in box
                    @test CartesianPoint{T}(box.x, box.y, box.z) in box
                    @test CartesianPoint{T}(-box.x, -box.y, -box.z) in box
                    @test !(CartesianPoint{T}(box.x, 0, 1.1*box.z) in box)
                    @test !(CartesianPoint{T}(box.y, box.y, box.y) in box)
                end
                @testset "Sphere" begin
                    sphere = Sphere(T(0.5)) # r from 0..0.5
                    @test CartesianPoint{T}(0, 0, 0) in sphere
                    @test CartesianPoint{T}(sphere.r, 0, 0) in sphere
                    @test CartesianPoint{T}(sphere.r/2, sphere.r/2, sphere.r/2) in sphere
                    @test !(CartesianPoint{T}(sphere.r, sphere.r, 0) in sphere)
                    @test !(CartesianPoint{T}(1.1*sphere.r, 0, 0) in sphere)
                    @test !(CartesianPoint{T}(0, 0.1*sphere.r, sphere.r) in sphere)
                end
                @testset "HexagonalPrism" begin
                    hexagon = HexagonalPrism(T(1.0), T(2.0), T(-1.0), T(1.0)) # outer radius from 1.0..2.0, z from -1.0..1.0
                    @test !(CartesianPoint{T}(0,0,0) in hexagon)
                    #in x: from sqrt(3)/2 to sqrt(3)
                    @test !(CartesianPoint{T}(0.5,0,0) in hexagon)
                    @test CartesianPoint{T}(sqrt(T(3))/2,0,0) in hexagon
                    @test CartesianPoint{T}(1,0,0) in hexagon
                    @test CartesianPoint{T}(sqrt(T(3)),0,0) in hexagon
                    @test !(CartesianPoint{T}(2,0,0) in hexagon)
                    #in y: from 1 to 2
                    @test !(CartesianPoint{T}(0,sqrt(T(3))/2,0) in hexagon)
                    @test CartesianPoint{T}(0,1,0) in hexagon
                    @test CartesianPoint{T}(0,sqrt(T(3)),0) in hexagon
                    @test CartesianPoint{T}(0,2,0) in hexagon
                    @test !(CartesianPoint{T}(0,2.01,0) in hexagon)
                    #in z: from -1 to 1
                    @test !(CartesianPoint{T}(0,1.5,-1.1) in hexagon)
                    @test CartesianPoint{T}(0,1.5,-1) in hexagon
                    @test CartesianPoint{T}(0,1.5,-0.5) in hexagon
                    @test CartesianPoint{T}(0,1.5,0.5) in hexagon
                    @test CartesianPoint{T}(0,1.5,1) in hexagon
                    @test !(CartesianPoint{T}(0,1.5,1.1) in hexagon)
                end
            end


            @testset "Transformations" begin
                # test transformations on a Tube
                tube = Tube(T(0), T(0.5), T(0), T(0), T(0), T(1)) # r from 0.0..0.5, z from 0.0..1.0

                @testset "Rotations" begin
                    RT = Float64
                    @testset "Rotation around X" begin
                        rotX = RotX{RT}(deg2rad(90))
                        rot_tube_x = RotMatrix3{T}(rotX) * tube # Tube with y from -1.0..0.0 and distance to y-axis of 0..0.5
                        @test CartesianPoint{T}(0, -1, 0) in rot_tube_x
                        @test CartesianPoint{T}(0, -0.5, 0) in rot_tube_x
                        @test CartesianPoint{T}(0, 0, 0) in rot_tube_x
                        @test !(CartesianPoint{T}(0, 0.5, 0) in rot_tube_x)
                        @test !(CartesianPoint{T}(0, 1, 0) in rot_tube_x)
                        @test CartesianPoint{T}(0.5, -0.5, 0) in rot_tube_x
                        @test CartesianPoint{T}(0, -0.5, 0.5) in rot_tube_x
                        @test CartesianPoint{T}(0.3, -0.5, 0.4) in rot_tube_x # should be exactly on the mantle
                        @test !(CartesianPoint{T}(0.5, -0.5, 0.5) in rot_tube_x)
                        @test !(CartesianPoint{T}(0, 0, 1) in rot_tube_x)
                    end
                    @testset "Single rotation around Y and Z" begin # first around Y, then around Z
                        rotZY = RotZY{RT}(deg2rad(-90), deg2rad(90))
                        rot_tube_zy = RotMatrix3{T}(rotZY) * tube # same Tube as rot_tube_x
                        @test CartesianPoint{T}(0, -1, 0) in rot_tube_zy
                        @test CartesianPoint{T}(0, -0.5, 0) in rot_tube_zy
                        @test CartesianPoint{T}(0, 0, 0) in rot_tube_zy
                        @test !(CartesianPoint{T}(0, 0.5, 0) in rot_tube_zy)
                        @test !(CartesianPoint{T}(0, 1, 0) in rot_tube_zy)
                        @test CartesianPoint{T}(0.5, -0.5, 0) in rot_tube_zy
                        @test CartesianPoint{T}(0, -0.5, 0.5) in rot_tube_zy
                        @test CartesianPoint{T}(0.3, -0.5, 0.4) in rot_tube_zy # should be exactly on the mantle
                        @test !(CartesianPoint{T}(0.5, -0.5, 0.5) in rot_tube_zy)
                        @test !(CartesianPoint{T}(0, 0, 1) in rot_tube_zy)
                    end
                    @testset "Repeated rotation around Y and Z" begin # first around Y, then around Z
                        rotZ = RotZ{RT}(deg2rad(-90))
                        rotY = RotY{RT}(deg2rad(90))
                        rot_tube_y =  RotMatrix3{T}(rotY) * tube
                        rot_tube_zy = RotMatrix3{T}(rotZ) * rot_tube_y # same Tube as before
                        @test CartesianPoint{T}(0, -1, 0) in rot_tube_zy
                        @test CartesianPoint{T}(0, -0.5, 0) in rot_tube_zy
                        @test CartesianPoint{T}(0, 0, 0) in rot_tube_zy
                        @test !(CartesianPoint{T}(0, 0.5, 0) in rot_tube_zy)
                        @test !(CartesianPoint{T}(0, 1, 0) in rot_tube_zy)
                        @test CartesianPoint{T}(0.5, -0.5, 0) in rot_tube_zy
                        @test CartesianPoint{T}(0, -0.5, 0.5) in rot_tube_zy
                        @test CartesianPoint{T}(0.3, -0.5, 0.4) in rot_tube_zy # should be exactly on the mantle
                        @test !(CartesianPoint{T}(0.5, -0.5, 0.5) in rot_tube_zy)
                        @test !(CartesianPoint{T}(0, 0, 1) in rot_tube_zy)
                    end
                    @testset "Single rotation around X, Y and Z" begin # first around X, then around Z, then around Y
                        rotYZX = RotYZX{RT}(deg2rad(90), deg2rad(90), deg2rad(90))
                        rot_tube_yzx = RotMatrix3{T}(rotYZX) * tube # Tube with z from -1.0..0.0 and distance to z-axis of 0..0.5
                        @test CartesianPoint{T}(0, 0, -1) in rot_tube_yzx
                        @test CartesianPoint{T}(0, 0, -0.5) in rot_tube_yzx
                        @test CartesianPoint{T}(0, 0, 0) in rot_tube_yzx
                        @test !(CartesianPoint{T}(0, 0, 0.5) in rot_tube_yzx)
                        @test !(CartesianPoint{T}(0, 0, 1) in rot_tube_yzx)
                        @test CartesianPoint{T}(0.5, 0, -1) in rot_tube_yzx
                        @test CartesianPoint{T}(0, -0.5, -1) in rot_tube_yzx
                        @test CartesianPoint{T}(0.4, 0.3, -0.5) in rot_tube_yzx # should be exactly on the mantle
                        @test !(CartesianPoint{T}(0.5, 0.5, -1) in rot_tube_yzx)
                        @test !(CartesianPoint{T}(0, 1, -0.5) in rot_tube_yzx)
                    end
                    @testset "Repeated rotation around X, Y and Z" begin # first around X, then around Z, then around Y
                        rotY = RotY{RT}(deg2rad(90))
                        rotZ = RotZ{RT}(deg2rad(90))
                        rotX = RotX{RT}(deg2rad(90))
                        rot_tube_x = RotMatrix3{T}(rotX) * tube
                        rot_tube_zx = RotMatrix3{T}(rotZ) * rot_tube_x
                        rot_tube_yzx = RotMatrix3{T}(rotY) * rot_tube_zx # same Tube as before
                        @test CartesianPoint{T}(0, 0, -1) in rot_tube_yzx
                        @test CartesianPoint{T}(0, 0, -0.5) in rot_tube_yzx
                        @test CartesianPoint{T}(0, 0, 0) in rot_tube_yzx
                        @test !(CartesianPoint{T}(0, 0, 0.5) in rot_tube_yzx)
                        @test !(CartesianPoint{T}(0, 0, 1) in rot_tube_yzx)
                        @test CartesianPoint{T}(0.5, 0, -1) in rot_tube_yzx
                        @test CartesianPoint{T}(0, -0.5, -1) in rot_tube_yzx
                        @test CartesianPoint{T}(0.4, 0.3, -0.5) in rot_tube_yzx # should be exactly on the mantle
                        @test !(CartesianPoint{T}(0.5, 0.5, -1) in rot_tube_yzx)
                        @test !(CartesianPoint{T}(0, 1, -0.5) in rot_tube_yzx)
                    end
                end

                @testset "Translations" begin
                    @testset "Single translation of a Tube along the axis" begin
                        translate1 = CartesianVector{T}(0, 0, -1)
                        translate_tube1 = tube + translate1 # same Tube as rot_tube_yzx
                        @test CartesianPoint{T}(0, 0, -1) in translate_tube1
                        @test CartesianPoint{T}(0, 0, -0.5) in translate_tube1
                        @test CartesianPoint{T}(0, 0, 0) in translate_tube1
                        @test !(CartesianPoint{T}(0, 0, 0.5) in translate_tube1)
                        @test !(CartesianPoint{T}(0, 0, 1) in translate_tube1)
                        @test CartesianPoint{T}(0.5, 0, -1) in translate_tube1
                        @test CartesianPoint{T}(0, -0.5, -1) in translate_tube1
                        @test !(CartesianPoint{T}(0.5, 0.5, -1) in translate_tube1)
                        @test !(CartesianPoint{T}(0, 1, -0.5) in translate_tube1)
                    end
                    @testset "Repeated translation of a Tube along the axis" begin
                        translate1 = CartesianVector{T}(0, 0, 1)
                        translate2 = CartesianVector{T}(0, 0, -2)
                        translate_tube = tube + translate1
                        translate_tube1 = translate_tube + translate2 # same Tube as before
                        @test CartesianPoint{T}(0, 0, -1) in translate_tube1
                        @test CartesianPoint{T}(0, 0, -0.5) in translate_tube1
                        @test CartesianPoint{T}(0, 0, 0) in translate_tube1
                        @test !(CartesianPoint{T}(0, 0, 0.5) in translate_tube1)
                        @test !(CartesianPoint{T}(0, 0, 1) in translate_tube1)
                        @test CartesianPoint{T}(0.5, 0, -1) in translate_tube1
                        @test CartesianPoint{T}(0, -0.5, -1) in translate_tube1
                        @test !(CartesianPoint{T}(0.5, 0.5, -1) in translate_tube1)
                        @test !(CartesianPoint{T}(0, 1, -0.5) in translate_tube1)
                    end
                    @testset "Single translation of a Tube perpendicular to the axis" begin
                        translate2 = CartesianVector{T}(1, 0, 0)
                        translate_tube2 = tube + translate2 # Tube with z from 0..1.0 and distance to z-axis at (1,0,0) of 0..0.5
                        @test !(CartesianPoint{T}(1, 0, -1) in translate_tube2)
                        @test !(CartesianPoint{T}(1, 0, -0.5) in translate_tube2)
                        @test CartesianPoint{T}(1, 0, 0) in translate_tube2
                        @test CartesianPoint{T}(1, 0, 0.5) in translate_tube2
                        @test CartesianPoint{T}(1, 0, 1) in translate_tube2
                        @test !(CartesianPoint{T}(0, 0, 0) in translate_tube2)
                        @test CartesianPoint{T}(0.5, 0, 0) in translate_tube2
                        @test CartesianPoint{T}(1.5, 0, 0) in translate_tube2
                        @test CartesianPoint{T}(1.3, -0.4, 0.5) in translate_tube2 # should be exactly on the mantle
                        @test !(CartesianPoint{T}(2, 0, 0) in translate_tube2)
                        @test !(CartesianPoint{T}(0.5, 0.5, 0) in translate_tube2)
                        @test CartesianPoint{T}(1, 0.5, 0) in translate_tube2
                    end
                    @testset "Repeated translation of a Tube perpendicular to the axis" begin
                        translate1 = CartesianVector{T}(0.5, 0.5, 0)
                        translate2 = CartesianVector{T}(0.5, -0.5, 0)
                        translate_tube1 = tube + translate1
                        translate_tube2 = translate_tube1 + translate2 # same Tube as before
                        @test !(CartesianPoint{T}(1, 0, -1) in translate_tube2)
                        @test !(CartesianPoint{T}(1, 0, -0.5) in translate_tube2)
                        @test CartesianPoint{T}(1, 0, 0) in translate_tube2
                        @test CartesianPoint{T}(1, 0, 0.5) in translate_tube2
                        @test CartesianPoint{T}(1, 0, 1) in translate_tube2
                        @test !(CartesianPoint{T}(0, 0, 0) in translate_tube2)
                        @test CartesianPoint{T}(0.5, 0, 0) in translate_tube2
                        @test CartesianPoint{T}(1.5, 0, 0) in translate_tube2
                        @test CartesianPoint{T}(1.3, -0.4, 0.5) in translate_tube2 # should be exactly on the mantle
                        @test !(CartesianPoint{T}(2, 0, 0) in translate_tube2)
                        @test !(CartesianPoint{T}(0.5, 0.5, 0) in translate_tube2)
                        @test CartesianPoint{T}(1, 0.5, 0) in translate_tube2
                    end
                end

                @testset "Scaling" begin
                    @testset "Single scaling of a Tube along the axis" begin
                        scale0 = SVector{3,T}(1.0, 1.0, 4.0)
                        scaled_tube = scale(tube, scale0) # Tube with z from 0..4.0 and distance to the z-axis of 0..0.5
                        @test !(CartesianPoint{T}(0, 0, -1) in scaled_tube)
                        @test CartesianPoint{T}(0, 0, 0) in scaled_tube
                        @test CartesianPoint{T}(0, 0, 1) in scaled_tube
                        @test CartesianPoint{T}(0, 0, 2) in scaled_tube
                        @test CartesianPoint{T}(0, 0, 3) in scaled_tube
                        @test CartesianPoint{T}(0, 0, 4) in scaled_tube
                        @test !(CartesianPoint{T}(0, 0, 5) in scaled_tube)
                        @test CartesianPoint{T}(0.5, 0, 4) in scaled_tube
                        @test CartesianPoint{T}(0, 0.5, 4) in scaled_tube
                        @test !(CartesianPoint{T}(0.5, 0.5, 4) in scaled_tube)
                    end
                    @testset "Repeated scaling of a Tube along the axis" begin
                        scale1 = SVector{3,T}(1.0, 1.0, 2.0)
                        scaled_tube1 = scale(tube, scale1)
                        scaled_tube2 = scale(scaled_tube1, scale1) # same Tube as before
                        @test !(CartesianPoint{T}(0, 0, -1) in scaled_tube2)
                        @test CartesianPoint{T}(0, 0, 0) in scaled_tube2
                        @test CartesianPoint{T}(0, 0, 1) in scaled_tube2
                        @test CartesianPoint{T}(0, 0, 2) in scaled_tube2
                        @test CartesianPoint{T}(0, 0, 3) in scaled_tube2
                        @test CartesianPoint{T}(0, 0, 4) in scaled_tube2
                        @test !(CartesianPoint{T}(0, 0, 5) in scaled_tube2)
                        @test CartesianPoint{T}(0.5, 0, 4) in scaled_tube2
                        @test CartesianPoint{T}(0, 0.5, 4) in scaled_tube2
                        @test !(CartesianPoint{T}(0.5, 0.5, 4) in scaled_tube2)
                    end
                    @testset "Single scaling of a Tube perpendicular to the axis" begin
                        scale3 = SVector{3,T}(2.0,0.5,1.0)
                        scaled_tube3 = scale(tube, scale3) # Tube with z from 0..1.0 and distance to the z-axis of 0..1.0 in x and 0..0.25 in y
                        @test CartesianPoint{T}(0, 0, 0) in scaled_tube3
                        @test CartesianPoint{T}(0.5, 0, 0) in scaled_tube3
                        @test CartesianPoint{T}(1, 0, 0) in scaled_tube3
                        @test !(CartesianPoint{T}(1.1, 0, 0) in scaled_tube3)
                        @test CartesianPoint{T}(0, 0, 0) in scaled_tube3
                        @test CartesianPoint{T}(0, 0.15, 0) in scaled_tube3
                        @test CartesianPoint{T}(0, 0.25, 0) in scaled_tube3
                        @test !(CartesianPoint{T}(0, 0.3, 0) in scaled_tube3)
                        @test CartesianPoint{T}(0.6, 0.2, 0) in scaled_tube3 #should be exactly on the mantle
                        @test !(CartesianPoint{T}(0.2, 0.6, 0) in scaled_tube3)
                    end
                    @testset "Repeated scaling of a Tube perpendicular to the axis" begin
                        scale4 = SVector{3,T}(2.0,1.0,1.0)
                        scale5 = SVector{3,T}(1.0,0.5,1.0)
                        scaled_tube4 = scale(tube, scale4)
                        scaled_tube5 = scale(scaled_tube4, scale5) # same Tube as before
                        @test CartesianPoint{T}(0, 0, 0) in scaled_tube5
                        @test CartesianPoint{T}(0.5, 0, 0) in scaled_tube5
                        @test CartesianPoint{T}(1, 0, 0) in scaled_tube5
                        @test !(CartesianPoint{T}(1.1, 0, 0) in scaled_tube5)
                        @test CartesianPoint{T}(0, 0, 0) in scaled_tube5
                        @test CartesianPoint{T}(0, 0.15, 0) in scaled_tube5
                        @test CartesianPoint{T}(0, 0.25, 0) in scaled_tube5
                        @test !(CartesianPoint{T}(0, 0.3, 0) in scaled_tube5)
                        @test CartesianPoint{T}(0.6, 0.2, 0) in scaled_tube5 #should be exactly on the mantle
                        @test !(CartesianPoint{T}(0.2, 0.6, 0) in scaled_tube5)
                    end
                end
            end

            @testset "Sets" begin
                @testset "Union" begin
                    @testset "Tube surrounded by Tube" begin
                        inner_tube = Tube(T(0.3),T(1.0)) #radius 0..0.3, height -0.5..0.5
                        outer_tube = Tube(T(0.7),T(1.5),T(0),T(0),T(-0.5),T(0.5)) #radius 0.7..1.0, height -0.5..0.5
                        union = inner_tube + outer_tube
                        @test CartesianPoint{T}(0.2,0.0,0.0) in union
                        @test CartesianPoint{T}(0.8,0.0,0.0) in union
                        @test !(CartesianPoint{T}(0.5,0.0,0.0) in union)
                        @test CartesianPoint{T}(0.2,0.0,0.1) in union
                        @test CartesianPoint{T}(0.8,0.0,0.1) in union
                        @test !(CartesianPoint{T}(0.5,0.0,0.1) in union)
                    end
                end
                @testset "Intersection" begin
                    @testset "Intersection between two Tubes" begin
                        tube1 = Tube(T(0),T(0.5),T(0),T(0),T(-1),T(0.5)) #radius 0..0.5, height -1.0..0.5
                        tube2 = Tube(T(0.3),T(1.0),T(0),T(0),T(-0.5),T(1.0)) #radius 0.3..1.0, height -0.5..1.0
                        intersection = tube1 & tube2 #should be tube with radius 0.3..0.5, height -0.5..0.5
                        @test !(CartesianPoint{T}(0.0,0.0,0.0) in intersection)
                        @test !(CartesianPoint{T}(0.2,0.0,0.0) in intersection)
                        @test CartesianPoint{T}(0.3,0.0,0.0) in intersection
                        @test CartesianPoint{T}(0.4,0.0,0.0) in intersection
                        @test CartesianPoint{T}(0.5,0.0,0.0) in intersection
                        @test !(CartesianPoint{T}(0.7,0.0,0.0) in intersection)
                        @test CartesianPoint{T}(0.3,0.0,0.5) in intersection
                        @test !(CartesianPoint{T}(0.4,0.0,-1.0) in intersection)
                        @test !(CartesianPoint{T}(0.5,0.0,1.0) in intersection)
                    end
                end
                @testset "Difference" begin
                    @testset "Difference between two Tubes" begin
                        inner_tube = Tube(T(0),T(0.5),T(0),T(0),T(-0.5),T(1)) # tube with radius 0..0.5, height -0.5..1.0
                        outer_tube = Tube(T(0.0),T(1.0),T(0),T(0),T(-1),T(1)) # tube with radius 0.0..1.0, height -1.0..1.0
                        difference = outer_tube - inner_tube # tube with radius 0.5..1.0, height -1.0..1.0
                        @test !(CartesianPoint{T}(0.0,0.0,0.0) in difference)
                        @test !(CartesianPoint{T}(0.5,0.0,0.0) in difference) #point on the surface should not be inside
                        @test CartesianPoint{T}(0.7,0.0,0.0) in difference
                        @test CartesianPoint{T}(1.0,0.0,0.0) in difference
                        @test !(CartesianPoint{T}(1.1,0.0,0.0) in difference)
                        @test !(CartesianPoint{T}(0.0,0.0,1.0) in difference)
                        #points of the subtracted volume should not be inside
                        @test !(CartesianPoint{T}(0.0,0.0,-0.5) in difference)
                        @test !(CartesianPoint{T}(0.5,0.0,-0.5) in difference)
                        @test CartesianPoint{T}(0.7,0.0,-0.5) in difference
                        @test CartesianPoint{T}(1.0,0.0,-0.5) in difference
                        #points on the surface of the initial volume should be inside if not subtracted
                        @test CartesianPoint{T}(0.0,0.0,-1.0) in difference
                        @test CartesianPoint{T}(0.5,0.0,-1.0) in difference
                        @test CartesianPoint{T}(0.7,0.0,-1.0) in difference
                        @test CartesianPoint{T}(1.0,0.0,-1.0) in difference
                    end
                end
            end
        end
    end
end
