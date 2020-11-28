using Test
using BenchmarkTools

using SolidStateDetectors
using IntervalSets
using Unitful
using LinearAlgebra
using Rotations
using StaticArrays

using SolidStateDetectors.ConstructiveSolidGeometry: 
    CartesianPoint, CartesianVector, CylindricalPoint,
    StretchedGeometry, RotatedGeometry, TranslatedGeometry,
    Tube, Cone, Box, Sphere, HexagonalPrism,
    CSGUnion, CSGIntersection, CSGDifference
    
code_warntype = true


T = Float64
p_car = CartesianPoint{T}(0, 1, 0)
@btime CylindricalPoint($p_car);


p_cyl = CylindricalPoint(p_car)
@btime CartesianPoint($p_cyl);

@testset "Points" begin
    T = Float64
    cp1 = CartesianPoint{T}(1, 2, 3)
    @test cp1 + cp1 == CartesianPoint{T}(2, 4, 6)
    @test cp1 - cp1 == CartesianPoint{T}(0, 0, 0)
    @test cp1 ⋅ cp1 == 14
    @test norm(cp1) == sqrt(14)
end


@testset "Volume Primitives" begin
    @testset "Tube" begin
        tube = Tube(0.2, 2.0) # radius 0.2, height 2
        @test CartesianPoint{T}(0, 0, 0) in tube
        @test CartesianPoint{T}(tube.r, 0, tube.z) in tube
        @test !(CartesianPoint{T}(tube.r, 0, 1.1*tube.z) in tube) 
        @test !(CartesianPoint{T}(1.1*tube.r, 0, tube.z) in tube) 
    end
    @testset "Cone" begin
        cone = Cone(0.5, 1.0, 2.0) # radius@top 0.5, radius@bottom 1, height 2
        @test CartesianPoint{T}(0, 0, 0) in cone
        @test CartesianPoint{T}(cone.rbot, 0, -cone.z) in cone
        @test !(CartesianPoint{T}(cone.rbot, 0, cone.z) in cone) 
        @test !(CartesianPoint{T}(0, 0, 1.1*cone.z) in cone) 
    end
    @testset "Box" begin
        box = Box(1.0,2.0,3.0) # x from -0.5..0.5, y from -1.0..1.0, z from -1.5..1.5
        @test CartesianPoint{T}(0, 0, 0) in box
        @test CartesianPoint{T}(box.x, box.y, box.z) in box
        @test CartesianPoint{T}(-box.x, -box.y, -box.z) in box
        @test !(CartesianPoint{T}(box.x, 0, 1.1*box.z) in box)
        @test !(CartesianPoint{T}(box.y, box.y, box.y) in box)
    end
    @testset "Sphere" begin
        sphere = Sphere(0.5) # radius 0.5
        @test CartesianPoint{T}(0, 0, 0) in sphere
        @test CartesianPoint{T}(sphere.r, 0, 0) in sphere
        @test CartesianPoint{T}(sphere.r/2, sphere.r/2, sphere.r/2) in sphere
        @test !(CartesianPoint{T}(sphere.r, sphere.r, 0) in sphere)
        @test !(CartesianPoint{T}(1.1*sphere.r, 0, 0) in sphere) 
        @test !(CartesianPoint{T}(0, 0.1*sphere.r, sphere.r) in sphere) 
    end
end


@testset "Stretching & Rotation & Translation of Primitives" begin 
    rot = RotX(deg2rad(90))
    t = CartesianVector{T}(2, 0, 0 )
    stretching = SVector(1.0, 3.0, 1.0)
    @testset "Tube" begin
        tube = Tube(0.2, 2.0) # r = 0.2, half_z = 1 -> h = 2 
        stretched_tube = StretchedGeometry(tube, stretching)
        rot_tube = RotatedGeometry(tube, RotMatrix(rot))
        translated_tube = TranslatedGeometry(tube, t)
        translated_rotated_tube = TranslatedGeometry(rot_tube, t)
        rot_stretched_tube =  RotatedGeometry(stretched_tube, RotMatrix(rot))
        translated_rot_stretched_tube = TranslatedGeometry(rot_stretched_tube, t)
        cp1 = CartesianPoint{T}(0, 0, 1)
        @test cp1 in tube
        @test cp1 in stretched_tube
        @test !(cp1 in rot_tube) 
        @test !(cp1 in rot_stretched_tube)
        @test !(cp1 in translated_tube) 
        @test !(cp1 in translated_rotated_tube) 
        @test !(cp1 in translated_rot_stretched_tube)
    end
    @testset "Cone" begin
        cone = Cone(0.5,1.0,2.0) # radius@top 0.5, radius@bottom 1, height 2
        cp1 = CartesianPoint{T}(0, 1, 0)
        stretched_cone = StretchedGeometry(cone, stretching)
        @test !(cp1 in cone)
        @test cp1 in stretched_cone
    end
end


rot = RotX(deg2rad(90))
t = CartesianVector{T}(1, 0, 0 )
tube = Tube(0, 0.1, 0, 2π, 0, 2.0) # r = 1, half_z = 1 -> h = 2 
stretching = SVector(1.0, 3.0, 1.0)
stretched_tube = StretchedGeometry(tube, stretching)
rot_tube = RotatedGeometry(tube, RotMatrix(rot))
translated_tube = TranslatedGeometry(tube, t)
translated_rotated_tube = TranslatedGeometry(rot_tube, t)
rot_stretched_tube =  RotatedGeometry(stretched_tube, RotMatrix(rot))
translated_rot_stretched_tube = TranslatedGeometry(rot_stretched_tube, t)
cp1 = CartesianPoint{T}(0, 0, 1)

if code_warntype
    @code_warntype cp1 in tube
    @code_warntype cp1 in stretched_tube
    @code_warntype cp1 in rot_tube
    @code_warntype cp1 in rot_stretched_tube
    @code_warntype cp1 in translated_tube
    @code_warntype cp1 in translated_rotated_tube
    @code_warntype cp1 in translated_rot_stretched_tube
end

@btime cp1 in translated_rot_stretched_tube


rot = RotX(deg2rad(90))
t = CartesianVector{T}(1, 0, 0 )
cone = Cone(0.1, 0.2, 0.1, 0.5, 0, 2π, 0, 2.0)
stretching = SVector(1.0, 3.0, 1.0)
stretched_cone = StretchedGeometry(cone, stretching)
rot_cone = RotatedGeometry(cone, RotMatrix(rot))
translated_cone = TranslatedGeometry(cone, t)
translated_rotated_cone = TranslatedGeometry(rot_cone, t)
rot_stretched_cone =  RotatedGeometry(stretched_cone, RotMatrix(rot))
translated_rot_stretched_cone = TranslatedGeometry(rot_stretched_cone, t)
cp1 = CartesianPoint{T}(0, 0, 1)

if code_warntype
    @code_warntype cp1 in cone
    @code_warntype cp1 in stretched_cone
    @code_warntype cp1 in rot_cone
    @code_warntype cp1 in rot_stretched_cone
    @code_warntype cp1 in translated_cone
    @code_warntype cp1 in translated_rotated_cone
    @code_warntype cp1 in translated_rot_stretched_cone
end

@btime cp1 in translated_rot_stretched_cone


rot = RotX(deg2rad(90))
t = CartesianVector{T}(1, 0, 0 )
box = Box(-1.0,1.0,-1.0,1.0,-1.0,1.0)
stretching = SVector(1.0, 3.0, 1.0)
stretched_box = StretchedGeometry(box, stretching)
rot_box = RotatedGeometry(box, RotMatrix(rot))
translated_box = TranslatedGeometry(box, t)
translated_rotated_box = TranslatedGeometry(rot_box, t)
rot_stretched_box =  RotatedGeometry(stretched_box, RotMatrix(rot))
translated_rot_stretched_box = TranslatedGeometry(rot_stretched_box, t)
cp1 = CartesianPoint{T}(0, 0, 1)

if code_warntype
    @code_warntype cp1 in box
    @code_warntype cp1 in stretched_box
    @code_warntype cp1 in rot_box
    @code_warntype cp1 in rot_stretched_box
    @code_warntype cp1 in translated_box
    @code_warntype cp1 in translated_rotated_box
    @code_warntype cp1 in translated_rot_stretched_box
end

@btime cp1 in translated_rot_stretched_box


rot = RotX(deg2rad(90))
t = CartesianVector{T}(1, 0, 0 )
sphere = Sphere(1.0)
stretching = SVector(1.0, 3.0, 1.0)
stretched_sphere = StretchedGeometry(sphere, stretching)
rot_sphere = RotatedGeometry(sphere, RotMatrix(rot))
translated_sphere = TranslatedGeometry(sphere, t)
translated_rotated_sphere = TranslatedGeometry(rot_sphere, t)
rot_stretched_sphere =  RotatedGeometry(stretched_sphere, RotMatrix(rot))
translated_rot_stretched_sphere = TranslatedGeometry(rot_stretched_sphere, t)
cp1 = CartesianPoint{T}(0, 0, 1)

if code_warntype
    @code_warntype cp1 in sphere
    @code_warntype cp1 in stretched_sphere
    @code_warntype cp1 in rot_sphere
    @code_warntype cp1 in rot_stretched_sphere
    @code_warntype cp1 in translated_sphere
    @code_warntype cp1 in translated_rotated_sphere
    @code_warntype cp1 in translated_rot_stretched_sphere
end

@btime cp1 in translated_rot_stretched_sphere


rot = RotX(deg2rad(90))
t = CartesianVector{T}(1, 0, 0 )
hexagon = HexagonalPrism(0.5,1.0,-1.0,1.0)
stretching = SVector(1.0, 3.0, 1.0)
stretched_hexagon = StretchedGeometry(hexagon, stretching)
rot_hexagon = RotatedGeometry(hexagon, RotMatrix(rot))
translated_hexagon = TranslatedGeometry(hexagon, t)
translated_rotated_hexagon = TranslatedGeometry(rot_hexagon, t)
rot_stretched_hexagon =  RotatedGeometry(stretched_hexagon, RotMatrix(rot))
translated_rot_stretched_hexagon = TranslatedGeometry(rot_stretched_hexagon, t)
cp1 = CartesianPoint{T}(0, 0, 1)

if code_warntype
    @code_warntype cp1 in hexagon
    @code_warntype cp1 in stretched_hexagon
    @code_warntype cp1 in rot_hexagon
    @code_warntype cp1 in rot_stretched_hexagon
    @code_warntype cp1 in translated_hexagon
    @code_warntype cp1 in translated_rotated_hexagon
    @code_warntype cp1 in translated_rot_stretched_hexagon
end

@btime cp1 in translated_rot_stretched_hexagon

using Plots
begin
    plts = []
    for p in [  
                tube, 
                stretched_tube, 
                rot_tube, 
                rot_stretched_tube, 
                translated_tube, 
                translated_rotated_tube ,
                translated_rot_stretched_tube,
            ]
        xs = T[]
        ys = T[]
        zs = T[]
        for x in -2:0.05:2
            for y in -2:0.05:2
                for z in -2:0.05:2
                    if CartesianPoint{T}(x, y, z) in p
                        push!(xs, x)
                        push!(ys, y)
                        push!(zs, z)
                    end
                end
            end
        end
        # push!(plts, plot(ys, zs, st=:scatter, ms = 4, markerstrokewidth = 0, xguide = "y", yguide = "z",
        #                     xlims = (-2,2), ylims = (-2,2)))
        push!(plts, plot3d(xs, ys, zs, st=:scatter, ms = 2, markerstrokewidth = 0, xguide = "x",
                            xlims = (-2,2), ylims = (-2,2), zlims = (-2,2)))
    end
    plot(plts..., size = (600,600))
end