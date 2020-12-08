using Test
using SolidStateDetectors
using IntervalSets
using LinearAlgebra
using JSON
using YAML

import SolidStateDetectors.ConstructiveSolidGeometry: Geometry
using SolidStateDetectors.ConstructiveSolidGeometry: 
    internal_unit_length, internal_unit_angle,
    Box, Cone, HexagonalPrism, Sphere, Tube

T = Float64

function Geometry(::Type{T}, filename::String, input_units::NamedTuple = (length = internal_unit_length, angle = internal_unit_angle)) where {T}
    @assert isfile(filename) "The given filename '$(filename)' does not lead to a valid file."
    dict = if endswith(filename, ".json")
        JSON.parsefile(filename)
    elseif endswith(filename, ".yaml")
        YAML.load_file(filename)
    else
        @error "Only JSON and YAML formats are supported at the moment."
    end
    Geometry(T, dict, input_units)
end


example_primitive_dir = joinpath(@__DIR__, "../../examples/example_primitive_files")
@testset "Test primitive read-in" begin
    @test typeof(Geometry(T, joinpath(example_primitive_dir, "Box.yaml"))) <: Box{T}
    @test typeof(Geometry(T, joinpath(example_primitive_dir, "Cone.yaml"))) <: Cone{T, <:Tuple}
    @test typeof(Geometry(T, joinpath(example_primitive_dir, "HexagonalPrism.yaml"))) <: HexagonalPrism{T}
    @test typeof(Geometry(T, joinpath(example_primitive_dir, "Sphere.yaml"))) <: Sphere{T}
    @test typeof(Geometry(T, joinpath(example_primitive_dir, "Tube.yaml"))) <: Cone{T, <:Union{T, AbstractInterval{T}}}
end