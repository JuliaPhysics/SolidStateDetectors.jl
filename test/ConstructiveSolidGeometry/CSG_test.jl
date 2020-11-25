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
    Tube,
    CSGUnion, CSGIntersection, CSGDifference


T = Float64
p_car = CartesianPoint{T}(0, 1, 0)

t1 = Tube(1, 0.4)
t2 = Tube(0.2, 1, 0, π, -0.3, 0.4)

p_car in t1
p_car in t2

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
        tube = Tube(0.2, 2.0) # radius 1, height 2
        @test CartesianPoint{T}(0, 0, 0) in tube
        @test CartesianPoint{T}(tube.r, 0, tube.z) in tube
        @test !(CartesianPoint{T}(tube.r, 0, 1.1*tube.z) in tube) 
        @test !(CartesianPoint{T}(1.1*tube.r, 0, tube.z) in tube) 
    end
end


@testset "Streching & Rotation & Translation of Primitives" begin 
    rot = RotX(deg2rad(90))
    t = CartesianVector{T}(2, 0, 0 )
    streching = SVector(1.0, 3.0, 1.0)

    tube = Tube(0.2, 2.0) # r = 1, half_z = 1 -> h = 2 

    streched_tube = StretchedGeometry(tube, streching)
    rot_tube = RotatedGeometry(tube, RotMatrix(rot))
    translated_tube = TranslatedGeometry(tube, t)
    translated_rotated_tube = TranslatedGeometry(rot_tube, t)
    rot_streched_tube =  RotatedGeometry(streched_tube, RotMatrix(rot))
    translated_rot_streched_tube = TranslatedGeometry(rot_streched_tube, t)

    cp1 = CartesianPoint{T}(0, 0, 1)
    
    @test (cp1 in tube)
    @test cp1 in streched_tube
    @test !(cp1 in rot_tube) 
    @test !(cp1 in rot_streched_tube)
    @test !(cp1 in translated_tube) 
    @test !(cp1 in translated_rotated_tube) 
    @test !(cp1 in translated_rot_streched_tube)
end
    
rot = RotX(deg2rad(90))
t = CartesianVector{T}(1, 0, 0 )
tube = Tube(0, 0.1, 0, 2π, 0, 2.0) # r = 1, half_z = 1 -> h = 2 
streching = SVector(1.0, 3.0, 1.0)
streched_tube = StretchedGeometry(tube, streching)
rot_tube = RotatedGeometry(tube, RotMatrix(rot))
translated_tube = TranslatedGeometry(tube, t)
translated_rotated_tube = TranslatedGeometry(rot_tube, t)
rot_streched_tube =  RotatedGeometry(streched_tube, RotMatrix(rot))
translated_rot_streched_tube = TranslatedGeometry(rot_streched_tube, t)
cp1 = CartesianPoint{T}(0, 0, 1)

@code_warntype cp1 in tube
@code_warntype cp1 in streched_tube
@code_warntype cp1 in rot_tube
@code_warntype cp1 in rot_streched_tube
@code_warntype cp1 in translated_tube
@code_warntype cp1 in translated_rotated_tube
@code_warntype cp1 in translated_rot_streched_tube

@btime cp1 in translated_rot_streched_tube

using Plots
begin
    plts = []
    for p in [  
                tube, 
                streched_tube, 
                rot_tube, 
                rot_streched_tube, 
                translated_tube, 
                translated_rotated_tube ,
                translated_rot_streched_tube,
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