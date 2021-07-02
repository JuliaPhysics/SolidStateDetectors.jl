using Test
using SolidStateDetectors
using IntervalSets
using LinearAlgebra
using JSON
using YAML
using StaticArrays

import SolidStateDetectors.ConstructiveSolidGeometry: Geometry
using SolidStateDetectors.ConstructiveSolidGeometry: 
    internal_unit_length, internal_unit_angle,
    Box, Cone, HexagonalPrism, Ellipsoid, Torus

T = Float64


example_primitive_dir = joinpath(@__DIR__, "../../examples/example_primitive_files")
@testset "Test primitive read-in" begin
    @test typeof(Geometry(T, joinpath(example_primitive_dir, "Box.yaml"))) <: Box{T}
    @test typeof(Geometry(T, joinpath(example_primitive_dir, "Cone.yaml"))) <: Cone{T, <:Any, <:Tuple}
    @test typeof(Geometry(T, joinpath(example_primitive_dir, "HexagonalPrism.yaml")).a) <: HexagonalPrism{T}
    @test typeof(Geometry(T, joinpath(example_primitive_dir, "Sphere.yaml")).a) <: Ellipsoid{T}
    @test typeof(Geometry(T, joinpath(example_primitive_dir, "Torus.yaml")).a) <: Torus{T}
    @test typeof(Geometry(T, joinpath(example_primitive_dir, "Tube.yaml"))) <: Cone{T, <:Any, <:Tuple}
end