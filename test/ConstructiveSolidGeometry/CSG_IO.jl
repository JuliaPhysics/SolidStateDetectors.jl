using Test
using SolidStateDetectors
using IntervalSets
using LinearAlgebra
using JSON
using YAML
using StaticArrays

import SolidStateDetectors.ConstructiveSolidGeometry as CSG
import SolidStateDetectors.ConstructiveSolidGeometry: Geometry
using SolidStateDetectors.ConstructiveSolidGeometry: 
    internal_unit_length, internal_unit_angle,
    Box, Cone, HexagonalPrism, Ellipsoid, Torus

T = Float64


example_primitive_dir = joinpath(@__DIR__, "../../examples/example_primitive_files")
@testset "Test primitive read-in" begin
    @testset "Box" begin
        box_widths = Geometry(T, joinpath(example_primitive_dir, "Box.yaml"))
        box_halfwidths = Geometry(T, joinpath(example_primitive_dir, "Box_halfwidths.yaml"))
        box_hXhYhZ = Geometry(T, joinpath(example_primitive_dir, "Box_hXhYhZ.yaml"))
        box_old = Geometry(T, joinpath(example_primitive_dir, "Box_old.yaml"))
        @test box_widths isa Box
        @test box_widths == box_halfwidths == box_hXhYhZ == box_old
    end

    @testset "Cone" begin
        cone = Geometry(T, joinpath(example_primitive_dir, "Cone.yaml"))
        cone_old = Geometry(T, joinpath(example_primitive_dir, "Cone_old.yaml"))
        @test cone isa Cone{T, <:Any, <:Tuple}
        @test cone == cone_old
    end

    @testset "Ellipsoid" begin
        ellipsoid_full_sphere = Geometry(T, joinpath(example_primitive_dir, "Ellipsoid_full_sphere.yaml"))
        @test ellipsoid_full_sphere isa CSG.FullSphere{T}
    end

    @test typeof(Geometry(T, joinpath(example_primitive_dir, "HexagonalPrism.yaml")).a) <: HexagonalPrism{T}
    @test typeof(Geometry(T, joinpath(example_primitive_dir, "Torus.yaml")).a) <: Torus{T}
    @test typeof(Geometry(T, joinpath(example_primitive_dir, "Tube.yaml"))) <: Cone{T, <:Any, <:Tuple}
end