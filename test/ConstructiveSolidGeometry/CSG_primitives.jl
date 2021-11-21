using Test
using SolidStateDetectors
using StaticArrays

import SolidStateDetectors.ConstructiveSolidGeometry as CSG
import SolidStateDetectors.ConstructiveSolidGeometry: Geometry, Dictionary

T = Float64

default_units = SolidStateDetectors.default_unit_tuple()
no_translations = (rotation = one(SMatrix{3, 3, T, 9}), translation = zero(CartesianVector{T}))

@testset "Test primitive parsing" begin
    @testset "Cone" begin
        for bot in (2.0, Dict("from" => 1.0, "to" => 2.0)),
            top in (2.0, 1.0, Dict("from" => 1.0, "to" => 2.0), Dict("from" => 2.0, "to" => 4.0)),
            φ in ((30,180),(0,360))
            
            dict = Dict("difference" => [
                Dict("cone"   => Dict(
                    "r"       => Dict("bottom" => bot, "top" => top),
                    "phi"     => Dict("from" => "$(φ[1])°", "to" => "$(φ[2])°"),
                    "h"       => 1.0))
                for i in 1:2
            ])
            
            c = Geometry(T, dict, default_units, no_translations)
            
            # Conversion from Geometry -> Dict and Dict -> Geometry should result in the same geometry
            output = Dictionary(c)
            name = collect(keys(output["difference"][1]))[1] # "tube" or "cone"
            if haskey(output["difference"][1][name], "phi")
                @test CSG._parse_value(T, output["difference"][1][name]["phi"]["to"], default_units.angle) == CSG._parse_value(T, "$(φ[2])°", default_units.angle)
            end
            # Internally, the primitive stores the offset of phi in a rotation, but it should not be part of the config file
            @test !haskey(output["difference"][1], "rotation")
            @test c.a.rotation ≈ Geometry(T, output, default_units, no_translations).a.rotation
            
            # No warnings or errors when decomposing the Cones into surfaces
            @test_nowarn CSG.surfaces(c.a)
            @test_nowarn CSG.surfaces(c.b)

            # Check if all Cones are saved the right way
            @test c.a isa CSG.Cone
            @test c.b isa CSG.Cone
        end
    end
    @testset "Torus" begin
        for r_tube in (2.0, Dict("from" => 1.0, "to" => 2.0)), 
            φ in ((30,180),(0,360)), 
            θmin in 0:45:90, 
            θmax in 135:45:270
            
            dict = Dict("difference" => [
                Dict( "torus" => Dict(
                    "r_torus" => 10.0,
                    "r_tube"  => r_tube,
                    "phi"     => Dict("from" => "$(φ[1])°", "to" => "$(φ[2])°"),
                    "theta"   => Dict("from" => "$(θmin)°", "to" => "$(θmax)°")))
                for i in 1:2
            ])
            t = Geometry(T, dict, default_units, no_translations)
            
            # Conversion from Geometry -> Dict and Dict -> Geometry should result in the same geometry
            output = Dictionary(t)
            if haskey(output["difference"][1]["torus"], "phi")
                @test CSG._parse_value(T, output["difference"][1]["torus"]["phi"]["to"], default_units.angle) == CSG._parse_value(T, "$(φ[2])°", default_units.angle)
            end
            # Internally, the primitive stores the offset of phi in a rotation, but it should not be part of the config file
            @test !haskey(output["difference"][1], "rotation") 
            @test t.a.rotation ≈ Geometry(T, output, default_units, no_translations).a.rotation
            
            # No warnings or errors when decomposing the Torus into surfaces
            @test_nowarn CSG.surfaces(t.a)
            @test_nowarn CSG.surfaces(t.b)

            # Check if all Torus are saved the right way
            @test t.a isa CSG.Torus
            @test t.b isa CSG.Torus
        end
    end
end