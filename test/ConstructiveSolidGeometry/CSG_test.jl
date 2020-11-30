using Test
using SolidStateDetectors
#using IntervalSets
#using Unitful
using LinearAlgebra
using Rotations
using StaticArrays

using SolidStateDetectors.ConstructiveSolidGeometry: 
    CartesianPoint, CartesianVector, CylindricalPoint,
    ScaledGeometry, RotatedGeometry, TranslatedGeometry,
    Tube, Cone, Box, Sphere, HexagonalPrism,
    CSGUnion, CSGIntersection, CSGDifference

T = Float64

@testset "Points" begin
    cp1 = CartesianPoint{T}(1, 2, 3)
    @test cp1 + cp1 == CartesianPoint{T}(2, 4, 6)
    @test cp1 - cp1 == CartesianPoint{T}(0, 0, 0)
    @test cp1 ⋅ cp1 == 14
    @test norm(cp1) == sqrt(14)
end

@testset "Volume Primitives" begin
    @testset "Tube" begin
        tube = Tube(0.2, 2.0) # r from 0..0.2, z from -1.0..1.0
        @test CartesianPoint{T}(0, 0, 0) in tube
        @test CartesianPoint{T}(tube.r, 0, tube.z) in tube
        @test !(CartesianPoint{T}(tube.r, 0, 1.1*tube.z) in tube) 
        @test !(CartesianPoint{T}(1.1*tube.r, 0, tube.z) in tube) 
    end
    @testset "Cone" begin
        cone = Cone(0.5, 1.0, 2.0) # radius@top from 0..0.5, radius@bottom from 0..1.0, z from -1.0..1.0
        @test CartesianPoint{T}(0, 0, 0) in cone
        @test CartesianPoint{T}(cone.rbot, 0, -cone.z) in cone
        @test !(CartesianPoint{T}(cone.rbot, 0, cone.z) in cone) 
        @test !(CartesianPoint{T}(0, 0, 1.1*cone.z) in cone) 
    end
    @testset "Box" begin
        box = Box(1.0, 2.0, 3.0) # x from -0.5..0.5, y from -1.0..1.0, z from -1.5..1.5
        @test CartesianPoint{T}(0, 0, 0) in box
        @test CartesianPoint{T}(box.x, box.y, box.z) in box
        @test CartesianPoint{T}(-box.x, -box.y, -box.z) in box
        @test !(CartesianPoint{T}(box.x, 0, 1.1*box.z) in box)
        @test !(CartesianPoint{T}(box.y, box.y, box.y) in box)
    end
    @testset "Sphere" begin
        sphere = Sphere(0.5) # r from 0..0.5
        @test CartesianPoint{T}(0, 0, 0) in sphere
        @test CartesianPoint{T}(sphere.r, 0, 0) in sphere
        @test CartesianPoint{T}(sphere.r/2, sphere.r/2, sphere.r/2) in sphere
        @test !(CartesianPoint{T}(sphere.r, sphere.r, 0) in sphere)
        @test !(CartesianPoint{T}(1.1*sphere.r, 0, 0) in sphere) 
        @test !(CartesianPoint{T}(0, 0.1*sphere.r, sphere.r) in sphere) 
    end
    @testset "HexagonalPrism" begin
        hexagon = HexagonalPrism(1.0, 2.0) # outer radius from 0..1.0, z from -1.0..1.0
        
    end
end


@testset "Transformations" begin 
    # test transformations on a Tube
    tube = Tube(0, 0.5, 0, 0, 0, 1) # r from 0.0..0.5, z from 0.0..1.0
    
    @testset "Rotations" begin
        @testset "Rotation around X" begin
            rotX = RotX{T}(deg2rad(90))
            rot_tube_x = RotatedGeometry(tube, RotMatrix(rotX)) # Tube with y from -1.0..0.0 and distance to y-axis of 0..0.5
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
        @testset "Rotation around Y and Z" begin # first around Y, then around Z
            rotZY = RotZY{T}(deg2rad(-90), deg2rad(90))
            rot_tube_zy = RotatedGeometry(tube, RotMatrix(rotZY)) # same Tube as rot_tube_x
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
        @testset "Rotations around X, Y and Z" begin # first around X, then around Z, then around Y
            rotYZX = RotYZX{T}(deg2rad(90), deg2rad(90), deg2rad(90))
            rot_tube_yzx = RotatedGeometry(tube, RotMatrix(rotYZX)) # Tube with z from -1.0..0.0 and distance to z-axis of 0..0.5
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
        @testset "Translation of a Tube along the axis" begin
            translate1 = CartesianVector{T}(0, 0, -1)
            translate_tube1 = TranslatedGeometry(tube, translate1) # same Tube as rot_tube_yzx
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
        @testset "Translation of a Tube perpendicular to the axis" begin
            translate2 = CartesianVector{T}(1, 0, 0)
            translate_tube2 = TranslatedGeometry(tube, translate2) # Tube with z from 0..1.0 and distance to z-axis at (1,0,0) of 0..0.5
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
        @testset "Scaling of a Tube along the axis" begin
            scale1 = SVector(1.0, 1.0, 4.0)
            scaled_tube1 = ScaledGeometry(tube, scale1) # Tube with z from 0..4.0 and distance to the z-axis of 0..0.5
            @test !(CartesianPoint{T}(0, 0, -1) in scaled_tube1)
            @test CartesianPoint{T}(0, 0, 0) in scaled_tube1
            @test CartesianPoint{T}(0, 0, 1) in scaled_tube1
            @test CartesianPoint{T}(0, 0, 2) in scaled_tube1
            @test CartesianPoint{T}(0, 0, 3) in scaled_tube1
            @test CartesianPoint{T}(0, 0, 4) in scaled_tube1
            @test !(CartesianPoint{T}(0, 0, 5) in scaled_tube1)
            @test CartesianPoint{T}(0.5, 0, 4) in scaled_tube1
            @test CartesianPoint{T}(0, 0.5, 4) in scaled_tube1
            @test !(CartesianPoint{T}(0.5, 0.5, 4) in scaled_tube1)
        end
        @testset "Scaling of a Tube perpendicular to the axis" begin
            scale2 = SVector(2.0,0.5,1.0)
            scaled_tube2 = ScaledGeometry(tube, scale2) # Tube with z from 0..1.0 and distance to the z-axis of 0..1.0 in x and 0..0.25 in y
            @test CartesianPoint{T}(0, 0, 0) in scaled_tube2
            @test CartesianPoint{T}(0.5, 0, 0) in scaled_tube2
            @test CartesianPoint{T}(1, 0, 0) in scaled_tube2
            @test !(CartesianPoint{T}(1.1, 0, 0) in scaled_tube2)
            @test CartesianPoint{T}(0, 0, 0) in scaled_tube2
            @test CartesianPoint{T}(0, 0.15, 0) in scaled_tube2
            @test CartesianPoint{T}(0, 0.25, 0) in scaled_tube2
            @test !(CartesianPoint{T}(0, 0.3, 0) in scaled_tube2)
            @test CartesianPoint{T}(0.6, 0.2, 0) in scaled_tube2 #should be exactly on the mantle
            @test !(CartesianPoint{T}(0.2, 0.6, 0) in scaled_tube2)
        end
    end
end

@testset "Sets" begin
    @testset "Union" begin
        @testset "Tube surrounded by Tube" begin
            inner_tube = Tube(0.3,1.0) #radius 0..0.3, height -0.5..0.5
            outer_tube = Tube(0.7,1.5,0,0,-0.5,0.5) #radius 0.7..1.0, height -0.5..0.5
            union = CSGUnion(inner_tube, outer_tube)
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
            tube1 = Tube(0,0.5,0,0,-1,0.5) #radius 0..0.5, height -1.0..0.5
            tube2 = Tube(0.3,1.0,0,0,-0.5,1.0) #radius 0.3..1.0, height -0.5..1.0
            intersection = CSGIntersection(tube1, tube2) #should be tube with radius 0.3..0.5, height -0.5..0.5
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
            inner_tube = Tube(0,0.5,0,0,-0.5,1) # tube with radius 0..0.5, height -0.5..1.0
            outer_tube = Tube(0.0,1.0,0,0,-1,1) # tube with radius 0.0..1.0, height -1.0..1.0
            difference = CSGDifference(outer_tube, inner_tube) # tube with radius 0.5..1.0, height -1.0..1.0
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
