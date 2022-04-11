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
    @testset "Ellipsoid" begin
        ellip1 = @inferred CSG.Ellipsoid(CSG.ClosedPrimitive,r=1f0,origin = zero(CartesianPoint{Float16}),rotation = one(SMatrix{3, 3, Float16, 9}))
        ellip2 = @inferred CSG.Ellipsoid{Float32}(r=1.0)
        @test ellip1 === ellip2
        
        dict = Dict("sphere"   => Dict(
                "r"       => 1.0))
        ellip3 = Geometry(T,dict,default_units,no_translations)
        output = Dictionary(ellip3)
        @test dict == output
        
        ellip_open = @inferred CSG.Ellipsoid(CSG.OpenPrimitive,r=1f0,origin = zero(CartesianPoint{Float16}),rotation = one(SMatrix{3, 3, Float16, 9}))
        @test in(CartesianPoint{Float32}(1,0,0),ellip1)
        @test !in(CartesianPoint{Float32}(1,0,0),ellip_open)
        
        ellip_closed_trafo = @inferred CSG.Ellipsoid(CSG.ClosedPrimitive,r=1.0, origin=1/sqrt(3)*CartesianPoint{Float32}(1,1,1),rotation=SMatrix{3}(0.5,sqrt(3)/2,0,-sqrt(3)/2,0.5,0,0,0,1))
        ellip_open_trafo = @inferred CSG.Ellipsoid(CSG.OpenPrimitive,r=1.0, origin=1/sqrt(3)*CartesianPoint{Float32}(1,1,1),rotation=SMatrix{3}(0.5,sqrt(3)/2,0,-sqrt(3)/2,0.5,0,0,0,1))
        @test !in(CartesianPoint{Float64}(-1e-8,0,0),ellip_closed_trafo)
        @test in(CartesianPoint{Float64}(1e-8,0,0),ellip_open_trafo)
    end
    @testset "Box" begin
        box1 = @inferred CSG.Box(CSG.ClosedPrimitive,hX=1f0, hY=2f0, hZ=1f0, origin = zero(CartesianPoint{Float16}),rotation = one(SMatrix{3, 3, Float16, 9}))
        box2 = @inferred CSG.Box{Float32}(hX=1.0, hY=2f0, hZ=1f0)
        @test box1 === box2
        
        dict = Dict("box"   => Dict(
                "hX"       => 1.0,
                "hY"       => 42.0,
                "hZ"       => 3.14))
        box3 = Geometry(T,dict,default_units,no_translations)
        output = Dictionary(box3)
        @test dict == output
        
        box_closed_trafo = @inferred CSG.Box(CSG.ClosedPrimitive, hX=2., origin=CartesianPoint{Float32}(1,1,1),rotation=SMatrix{3}(0,0,1,0,1,0,-1,0,0))
        box_open_trafo = @inferred CSG.Box(CSG.OpenPrimitive, hX=2., origin=CartesianPoint{Float32}(1,1,1),rotation=SMatrix{3}(0,0,1,0,1,0,-1,0,0))
        tol=1e-8
        @test in(CartesianPoint{Float64}(1,1,3),box_closed_trafo)
        @test !in(CartesianPoint{Float64}(1,1,3),box_open_trafo)
        @test !in(CartesianPoint{Float64}(1,1,3+tol),box_closed_trafo)
        @test in(CartesianPoint{Float64}(1,1,3-tol),box_open_trafo)
    end
end