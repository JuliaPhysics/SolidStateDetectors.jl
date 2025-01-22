using Test
using SolidStateDetectors

import SolidStateDetectors.ConstructiveSolidGeometry as CSG
import SolidStateDetectors.ConstructiveSolidGeometry: Geometry

T = Float64


example_primitive_dir = joinpath(@__DIR__, "../../examples/example_primitive_files")
@testset "Test primitive read-in" begin
    @testset "Box" begin
        box_widths = Geometry(T, joinpath(example_primitive_dir, "Box.yaml"))
        box_halfwidths = Geometry(T, joinpath(example_primitive_dir, "Box_halfwidths.yaml"))
        box_hXhYhZ = Geometry(T, joinpath(example_primitive_dir, "Box_hXhYhZ.yaml"))
        @test box_widths isa CSG.Box
        @test box_widths == box_halfwidths == box_hXhYhZ
    end

    @testset "Cone" begin
        cone = Geometry(T, joinpath(example_primitive_dir, "Cone.yaml"))
        @test cone isa CSG.Cone{T, <:Any, <:Tuple}

        cone_tube = Geometry(T, joinpath(example_primitive_dir, "Cone_tube.yaml"))
        @test cone_tube isa CSG.VaryingTube{T}
    end

    @testset "Ellipsoid" begin
        ellipsoid_full_sphere = Geometry(T, joinpath(example_primitive_dir, "Ellipsoid_full_sphere.yaml"))
        @test ellipsoid_full_sphere isa CSG.FullSphere{T}
    end

    @testset "Polycone" begin
        polycone = Geometry(T, joinpath(example_primitive_dir, "Polycone.yaml"))
        @test polycone isa CSG.Polycone{T}
    end

    @testset "RegularPrism" begin
        hexagon = Geometry(T, joinpath(example_primitive_dir, "RegularPrism_hexagon.yaml"))
        @test hexagon isa CSG.HexagonalPrism{T}
    end

    @testset "Torus" begin
        torus = Geometry(T, joinpath(example_primitive_dir, "Torus.yaml"))
        @test torus isa CSG.Torus{T}
    end
end