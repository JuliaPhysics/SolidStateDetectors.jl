using Test
using SolidStateDetectors
using StaticArrays
using LinearAlgebra
using Rotations
using Unitful

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
    
        # Test different constructor versions and dictionary construction
        cone1 = @inferred CSG.Cone(CSG.ClosedPrimitive,r=1f0, φ = π, hZ = 1f0, origin = zero(CartesianPoint{Float16}),rotation = one(SMatrix{3, 3, Float16, 9}))
        cone2 = @inferred CSG.Cone{Float32}(r=1.0, φ = π)
        @test cone1 === cone2
        dict = Dict("tube"   => Dict(
                "r"       => 1.0,
                "h"       => 1.0))
        cone4 = Geometry(T,dict,default_units,no_translations)
        output = Dictionary(cone4)
        @test dict == output
    
        ## Test in method for several cones
        # Non transformed geometries
        @test !in(CartesianPoint{Float32}(0,-1,0),cone1)
        @test in(CartesianPoint{Float64}(0,-1,0),cone4)
        cone_open = @inferred CSG.Cone(CSG.OpenPrimitive, r=1f0,hZ=1f0, origin = zero(CartesianPoint{Float32}),rotation = one(SMatrix{3, 3, Float32, 9}))
        @test !in(CartesianPoint{Float32}(1,0,0),cone_open)
        # Boundary behaviour of closed and open primitives for transformed 
        # geometries (origin, rotation)
        cone_closed_trafo = @inferred CSG.Cone(CSG.ClosedPrimitive,r=1.0, hZ=3.0, origin=CartesianPoint{Float32}(0,0,1),rotation=SMatrix{3}(0,0,-1,0,1,0,1,0,0))#90 degree rotation arond y-axis
        cone_open_trafo = @inferred CSG.Cone(CSG.OpenPrimitive,r=1.0, hZ=3.0, origin=CartesianPoint{Float32}(0,0,1),rotation=SMatrix{3}(0,0,-1,0,1,0,1,0,0))
        @test CSG.isClosedPrimitive(cone_closed_trafo)
        @test !CSG.isClosedPrimitive(cone_open_trafo)
        @test cone_closed_trafo == CSG.switchClosedOpen(cone_open_trafo)
        @test cone_open_trafo == CSG.switchClosedOpen(cone_closed_trafo)
        @test in(CartesianPoint{Float64}(3,0,1),cone_closed_trafo)
        @test !in(CartesianPoint{Float64}(3,0,1),cone_open_trafo)
        tol=1e-8
        @test !in(CartesianPoint{Float64}(3+tol,0,1),cone_closed_trafo)
        @test in(CartesianPoint{Float64}(3-tol,0,1),cone_open_trafo)
        #Cones with r being any valid combination
        cone_empty = @inferred CSG.Cone{Float64}(r=((1,2),(3,5)), hZ=2, rotation = SMatrix{3}(0,0,-1,0,1,0,1,0,0))
        @test in(CartesianPoint{Float64}(1,0,2.5+tol), cone_empty)
        @test !in(CartesianPoint{Float64}(1,0,2.5-tol), cone_empty)
        @inferred CSG.Cone{Float64}(r=(nothing,(3,5)))
        @inferred CSG.Cone{Float64}(r=((3,5),nothing))
        @inferred CSG.Cone{Float64}(r=(1,2))
        @inferred CSG.Cone{Float64}(r=(nothing,(3,5)))
        @inferred CSG.Cone{Float64}(r=((1,),(2,)))
        @inferred CSG.Cone(r=(nothing,(3,5)))
        @inferred CSG.Cone(r=((3,5),nothing))
        @inferred CSG.Cone(r=(1,2))
        @inferred CSG.Cone(r=(nothing,(3,5)))
        @inferred CSG.Cone(r=((1,),(2,)))
        #Test where phi is a Tuple
        tuple_cone = @inferred CSG.Cone{Float64}(φ=(π/4,3*π/4), rotation=SMatrix{3}(0,0,-1,0,1,0,1,0,0))
        rot_cone = @inferred CSG.Cone{Float64}(φ=π/2, rotation=SMatrix{3}(0,0,-1,0,1,0,1,0,0) * SMatrix{3}(cos(π/4),sin(π/4),0,-sin(π/4),cos(π/4),0,0,0,1))
        @test tuple_cone ==rot_cone
        
        # PartialCylinder (OpenPrimitive)
        pc_open = CSG.Cone(CSG.OpenPrimitive, r = convert(T, 1.0), φ = convert(T, π/2), hZ = convert(T, 1.0))
        @test CSG._in(CartesianPoint{T}(0.5, 0.2, 0.0), pc_open)
        @test !CSG._in(CartesianPoint{T}(1.1, 0.0, 0.0), pc_open)
        @test !CSG._in(CartesianPoint{T}(-0.5, -0.2, 0.0), pc_open)
        @test !CSG._in(CartesianPoint{T}(0.5, 0.0, 1.0), pc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 0.0, 0.0), pc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 0.0, -1.0), pc_open)
        @test !CSG._in(CartesianPoint{T}(1.0, 0.0, -1.0), pc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 1.0, -1.0), pc_open)
        @test !CSG._in(CartesianPoint{T}(-1.0, 0.0, -1.0), pc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 0.0, 1.0), pc_open)
        
        # PartialCylinder (ClosedPrimitive) 
        pc_closed = CSG.Cone(CSG.ClosedPrimitive, r = convert(T, 1.0), φ = convert(T, π/2), hZ = convert(T, 1.0))
        @test Geometry(T, Dictionary(pc_closed), default_units, no_translations) == pc_closed
        @test CSG._in(CartesianPoint{T}(0.5, 0.2, 0.0), pc_closed)
        @test !CSG._in(CartesianPoint{T}(1.1, 0.0, 0.0), pc_closed)
        @test !CSG._in(CartesianPoint{T}(-0.5, -0.2, 0.0), pc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 0.0, 0.0), pc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 0.0, -1.0), pc_closed)
        @test CSG._in(CartesianPoint{T}(1.0, 0.0, -1.0), pc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 1.0, -1.0), pc_closed)
        @test !CSG._in(CartesianPoint{T}(-1.0, 0.0, -1.0), pc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 0.0, 1.0), pc_closed)
        
        # VaryingCylinder (ClosedPrimitive)
        vc_closed = CSG.Cone{T}(r = ((convert(T,1.0),), (convert(T,3.0),)), φ = nothing, hZ = convert(T,2.0))
        @test Geometry(T, Dictionary(vc_closed), default_units, no_translations) == vc_closed
        @test CSG._in(CartesianPoint{T}(1.9, 0.0, 0.0), vc_closed)
        @test !CSG._in(CartesianPoint{T}(2.1, 0.0, 0.0), vc_closed)
        @test CSG._in(CartesianPoint{T}(1.0, 0.0, -2.0), vc_closed)
        @test CSG._in(CartesianPoint{T}(3.0, 0.0, 2.0), vc_closed)
        
        # VaryingCylinder (OpenPrimitive)
        vc_open = CSG.Cone(CSG.OpenPrimitive, r = ((convert(T,1.0),), (convert(T,3.0),)), φ = nothing, hZ = convert(T,2.0))
        @test !CSG._in(CartesianPoint{T}(2.0, 0.0, 0.0), vc_open)
        @test CSG._in(CartesianPoint{T}(1.9, 0.0, 0.0), vc_open)
        @test !CSG._in(CartesianPoint{T}(1.0, 0.0, -2.0), vc_open)
        @test !CSG._in(CartesianPoint{T}(3.0, 0.0, 2.0), vc_open)
        
        # PartialVaryingCylinder (ClosedPrimitive)
        pvc_closed = CSG.Cone{T}(CSG.ClosedPrimitive, r = ((convert(T,1.0),), (convert(T,2.0),)), φ = convert(T, π), hZ = convert(T, 1.0))
        @test Geometry(T, Dictionary(pvc_closed), default_units, no_translations) == pvc_closed
        @test CSG._in(CartesianPoint{T}(1.4, 0.1, 0.0), pvc_closed)
        @test !CSG._in(CartesianPoint{T}(1.6, 0.0, 0.0), pvc_closed)
        @test !CSG._in(CartesianPoint{T}(1.0, -0.1, 0.0), pvc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 0.0, -1.0), pvc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 1.0, -1.0), pvc_closed)
        @test CSG._in(CartesianPoint{T}(1.0, 0.0, -1.0), pvc_closed)
        @test CSG._in(CartesianPoint{T}(-1.0, 0.0, -1.0), pvc_closed)
        @test !CSG._in(CartesianPoint{T}(0.0, -1.0, -1.0), pvc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 0.0, 1.0), pvc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 0.0, 0.0), pvc_closed)
        @test CSG._in(CartesianPoint{T}(1.5, 0.0, 1.0), pvc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 2.0, 1.0), pvc_closed)
        @test CSG._in(CartesianPoint{T}(2.0, 0.0, 1.0), pvc_closed)
        @test CSG._in(CartesianPoint{T}(-2.0, 0.0, 1.0), pvc_closed)
        @test !CSG._in(CartesianPoint{T}(0.0, -2.0, 1.0), pvc_closed)
        
        # PartialVaryingCylinder (OpenPrimitive)
        pvc_open = CSG.Cone(CSG.OpenPrimitive, r = ((convert(T,1.0),), (convert(T,2.0),)), φ = convert(T, π), hZ = convert(T, 1.0))
        @test !CSG._in(CartesianPoint{T}(1.5, 0.0, 0.0), pvc_open)
        @test CSG._in(CartesianPoint{T}(1.4, 0.1, 0.0), pvc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 0.0, -1.0), pvc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 1.0, -1.0), pvc_open)
        @test !CSG._in(CartesianPoint{T}(1.0, 0.0, -1.0), pvc_open)
        @test !CSG._in(CartesianPoint{T}(-1.0, 0.0, -1.0), pvc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, -1.0, -1.0), pvc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 0.0, 1.0), pvc_open)
        @test !CSG._in(CartesianPoint{T}(1.5, 0.0, 1.0), pvc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 2.0, 1.0), pvc_open)
        @test !CSG._in(CartesianPoint{T}(2.0, 0.0, 1.0), pvc_open)
        @test !CSG._in(CartesianPoint{T}(-2.0, 0.0, 1.0), pvc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, -2.0, 1.0), pvc_open)
        
        # PartialVaryingTube (ClosedPrimitive)
        pvt_closed = CSG.Cone(CSG.ClosedPrimitive, r = ((convert(T,1.0), convert(T,2.0)), (convert(T,1.5), convert(T,2.5))), φ = convert(T, π/2), hZ = convert(T, 1.0))
        @test Geometry(T, Dictionary(pvt_closed), default_units, no_translations) == pvt_closed
        @test CSG._in(CartesianPoint{T}(1.6, 0.2, 0.0), pvt_closed)
        @test !CSG._in(CartesianPoint{T}(1.1, 0.0, 0.0), pvt_closed)
        @test !CSG._in(CartesianPoint{T}(2.6, 0.0, 0.0), pvt_closed)
        @test !CSG._in(CartesianPoint{T}(1.6, -0.2, 0.0), pvt_closed)
        @test CSG._in(CartesianPoint{T}(1.0, 0.0, -1.0), pvt_closed)
        @test CSG._in(CartesianPoint{T}(2.0, 0.0, -1.0), pvt_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 1.0, -1.0), pvt_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 2.0, -1.0), pvt_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 2.5, 1.0), pvt_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 1.5, 1.0), pvt_closed)
        @test CSG._in(CartesianPoint{T}(1.5, 0.0, 1.0), pvt_closed)
        @test CSG._in(CartesianPoint{T}(2.5, 0.0, 1.0), pvt_closed)
        @test !CSG._in(CartesianPoint{T}(-2.5, 0.0, 1.0), pvt_closed)
        
        # PartialVaryingTube (OpenPrimitive)
        pvt_open = CSG.Cone(CSG.OpenPrimitive, r = ((convert(T,1.0), convert(T,2.0)), (convert(T,1.5), convert(T,2.5))), φ = convert(T, π/2), hZ = convert(T, 1.0))
        @test CSG._in(CartesianPoint{T}(1.6, 0.2, 0.0), pvt_open)
        @test !CSG._in(CartesianPoint{T}(1.1, 0.0, 0.0), pvt_open)
        @test !CSG._in(CartesianPoint{T}(2.6, 0.0, 0.0), pvt_open)
        @test !CSG._in(CartesianPoint{T}(1.6, -0.2, 0.0), pvt_open)
        @test !CSG._in(CartesianPoint{T}(1.0, 0.0, -1.0), pvt_open)
        @test !CSG._in(CartesianPoint{T}(2.0, 0.0, -1.0), pvt_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 1.0, -1.0), pvt_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 2.0, -1.0), pvt_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 2.5, 1.0), pvt_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 1.5, 1.0), pvt_open)
        @test !CSG._in(CartesianPoint{T}(1.5, 0.0, 1.0), pvt_open)
        @test !CSG._in(CartesianPoint{T}(2.5, 0.0, 1.0), pvt_open)
        @test !CSG._in(CartesianPoint{T}(-2.5, 0.0, 1.0), pvt_open)
        
        # Upward cones
        r_up_closed = ((nothing, convert(T,0.1)), nothing)
        r_up_open   = ((nothing, convert(T,0.1)), nothing)
        φ_up = nothing
        h_up = convert(T, 1.0)
        origin = zero(CartesianPoint{Float64}())
        rotation = SMatrix{3,3,T}(I)
        
        # Closed UpwardCone
        upward_closed = CSG.UpwardCone{T, CSG.ClosedPrimitive}(r_up_closed, φ_up, h_up, origin, rotation)
        @test Geometry(T, Dictionary(upward_closed), default_units, no_translations) == upward_closed
        @test CSG._in(CartesianPoint{T}(0.0, 0.0, 0.5), upward_closed)
        @test !CSG._in(CartesianPoint{T}(0.0, 0.0, 1.1), upward_closed)
        @test !CSG._in(CartesianPoint{T}(2.0, 0.0, 0.5), upward_closed)
        @test length(CSG.surfaces(upward_closed)) > 0
        
        # Open UpwardCone
        upward_open = CSG.UpwardCone{T, CSG.OpenPrimitive}(r_up_open, φ_up, h_up, origin, rotation)
        @test CSG._in(CartesianPoint{T}(0.0, 0.0, 0.5), upward_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 0.0, 1.0), upward_open)
        @test !CSG._in(CartesianPoint{T}(2.0, 0.0, 0.5), upward_open)
        @test length(CSG.surfaces(upward_open)) > 0
        
        r_bot = (nothing, 2.0)
        r_cone = (r_bot, nothing)
        
        # PartialUpwardCone 
        puc_closed = CSG.PartialUpwardCone{T, CSG.ClosedPrimitive}(r_cone, T(π), 1.0, origin, rotation)
        @test Geometry(T, Dictionary(puc_closed), default_units, no_translations) == puc_closed
        @test CSG._in(CartesianPoint{T}(0.5, 0.5, -0.5), puc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 0.0, -1.0), puc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 0.0, 0.0), puc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 0.0, 1.0), puc_closed)
        @test CSG._in(CartesianPoint{T}(2.0, 0.0, -1.0), puc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 2.0, -1.0), puc_closed)
        @test CSG._in(CartesianPoint{T}(-2.0, 0.0, -1.0), puc_closed)
        @test !CSG._in(CartesianPoint{T}(0.0, -2.0, -1.0), puc_closed)
        @test !CSG._in(CartesianPoint{T}(0.5, -0.6, -2), puc_closed)
        
        puc_open   = CSG.PartialUpwardCone{T, CSG.OpenPrimitive}(r_cone, T(π), 1.0, origin, rotation)
        @test CSG._in(CartesianPoint{T}(0.5, 0.5, -0.5), puc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 0.0, -1.0), puc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 0.0, 0.0), puc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 0.0, 1.0), puc_open)
        @test !CSG._in(CartesianPoint{T}(2.0, 0.0, -1.0), puc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 2.0, -1.0), puc_open)
        @test !CSG._in(CartesianPoint{T}(-2.0, 0.0, -1.0), puc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, -2.0, -1.0), puc_open)
        @test !CSG._in(CartesianPoint{T}(0.5, -0.6, -2), puc_open)
        
        surfs = CSG.surfaces(puc_open)
        @test length(surfs) > 0
        @test all(x -> x !== nothing, surfs)
        
        # Radii must match type: Tuple{Nothing, Tuple{Nothing,T}}
        r_top = (nothing, 2.0)
        r_cone = (nothing, r_top)
                
        # DownwardCone ClosedPrimitive
        dc_closed = CSG.DownwardCone{T, CSG.ClosedPrimitive}(r_cone, nothing, 1.0, origin, rotation)
        @test Geometry(T, Dictionary(dc_closed), default_units, no_translations) == dc_closed
        @test CSG._in(CartesianPoint{T}(0.2, 0.2, -0.5), dc_closed)
        @test CSG._in(CartesianPoint{T}(1, 0, 0), dc_closed)
        @test !CSG._in(CartesianPoint{T}(1.5, 0.0, -0.5), dc_closed)
        @test CSG._in(CartesianPoint{T}(0, 0, -1), dc_closed)
        
        # DownwardCone OpenPrimitive
        dc_open = CSG.DownwardCone{T, CSG.OpenPrimitive}(r_cone, nothing, 1.0, origin, rotation)
        @test CSG._in(CartesianPoint{T}(0.2, 0.2, -0.5), dc_open)
        @test !CSG._in(CartesianPoint{T}(1, 0, 0), dc_open)
        @test !CSG._in(CartesianPoint{T}(1.5, 0.0, -0.5), dc_open)
        @test !CSG._in(CartesianPoint{T}(0, 0, -1), dc_open)
        
        # PartialDownwardCone ClosedPrimitive
        pdc_closed = CSG.PartialDownwardCone{T, CSG.ClosedPrimitive}(r_cone, T(π), 1.0, origin, rotation)
        @test Geometry(T, Dictionary(pdc_closed), default_units, no_translations) == pdc_closed
        @test CSG._in(CartesianPoint{T}(0.5, 0.0, -0.5), pdc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 0.0, -1.0), pdc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 0.0, 0.0), pdc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 0.0, 1.0), pdc_closed)
        @test CSG._in(CartesianPoint{T}(2.0, 0.0, 1.0), pdc_closed)
        @test CSG._in(CartesianPoint{T}(0.0, 2.0, 1.0), pdc_closed)
        @test CSG._in(CartesianPoint{T}(-2.0, 0.0, 1.0), pdc_closed)
        @test !CSG._in(CartesianPoint{T}(0.0, -2.0, 1.0), pdc_closed)
        @test !CSG._in(CartesianPoint{T}(0.5, -0.6, -0.5), pdc_closed)
        
        # PartialDownwardCone OpenPrimitive
        pdc_open = CSG.PartialDownwardCone{T, CSG.OpenPrimitive}(r_cone, T(π), 1.0, origin, rotation)
        @test CSG._in(CartesianPoint{T}(0.2, 0.2, -0.5), pdc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 0.0, -1.0), pdc_open)
        @test !CSG._in(CartesianPoint{T}(1.5, 0.0, -0.5), pdc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 0.0, 0.0), pdc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 0.0, 1.0), pdc_open)
        @test !CSG._in(CartesianPoint{T}(2.0, 0.0, 1.0), pdc_open)
        @test !CSG._in(CartesianPoint{T}(0.0, 2.0, 1.0), pdc_open)
        @test !CSG._in(CartesianPoint{T}(-2.0, 0.0, 1.0), pdc_open)
        
        # Surfaces tests
        surfs_dc_closed = CSG.surfaces(dc_closed)
        @test length(surfs_dc_closed) == 2
        @test surfs_dc_closed[1] isa CSG.EllipticalSurface
        @test surfs_dc_closed[2] isa CSG.FullConeMantle
        
        surfs_pdc_closed = CSG.surfaces(pdc_closed)
        @test length(surfs_pdc_closed) == 4
        @test surfs_pdc_closed[1] isa CSG.EllipticalSurface
        @test surfs_pdc_closed[2] isa CSG.PartialConeMantle
        @test surfs_pdc_closed[3] isa CSG.Triangle
        @test surfs_pdc_closed[4] isa CSG.Triangle
        
        # TopClosedTube
        r_inner = (0.3, 0.5)
        r_top = 1.0
        r_cone = (r_inner, r_top)
        hZ = 1.0
        tube_closed = CSG.TopClosedTube{Float64, CSG.ClosedPrimitive}(r_cone, nothing, hZ, origin, rotation)
        @test CSG._in(CartesianPoint(0.5, 0.0, -1), tube_closed)    # surface limit
        @test CSG._in(CartesianPoint(1.0, 0.0, 1.0), tube_closed)   # on outer boundary
        @test !CSG._in(CartesianPoint(0.0, 0.0, -1.5), tube_closed) # below bottom
        
        tube_open = CSG.TopClosedTube{Float64, CSG.OpenPrimitive}(r_cone, nothing, hZ, origin, rotation)
        @test !CSG._in(CartesianPoint(0.5, 0.0, -1), tube_open)
        @test !CSG._in(CartesianPoint(1.0, 0.0, 1.0), tube_open)
        @test !CSG._in(CartesianPoint(0.0, 0.0, -1.5), tube_open)
        
        surfs_closed = CSG.surfaces(tube_closed)
        @test length(surfs_closed) == 3
        @test surfs_closed[1] isa CSG.EllipticalSurface
        @test surfs_closed[2] isa CSG.FullConeMantle
        @test surfs_closed[3] isa CSG.FullConeMantle
        
        ε = 1e-6
        pt = surfs_closed[1].origin + CartesianVector(r_inner[2], 0.0, 0.0)
        n = CSG.normal(surfs_closed[1], pt)
        pt_plus  = pt + ε*n
        pt_minus = pt - ε*n
        @test CSG._in(pt_plus, tube_closed)
        @test !CSG._in(pt_minus, tube_closed)
        
        surfs_open = CSG.surfaces(tube_open)
        @test length(surfs_open) == 3
        @test surfs_open[1] isa CSG.EllipticalSurface
        @test surfs_open[2] isa CSG.FullConeMantle
        @test surfs_open[3] isa CSG.FullConeMantle
        
        phi_partial = π/2
        
        # PartialTopClosedTube (ClosedPrimitive)
        p_tube_closed = CSG.PartialTopClosedTube{Float64, CSG.ClosedPrimitive}(r_cone, phi_partial, hZ, origin, rotation)
        @test CSG._in(CartesianPoint(1.0, 0.0, 1.0), p_tube_closed)
        @test CSG._in(CartesianPoint(1.0, 0.0, 1.0), p_tube_closed)
        @test !CSG._in(CartesianPoint(-1.0, 0.0, 1.0), p_tube_closed)
        @test CSG._in(CartesianPoint(0.5, 0.0, -1.0), p_tube_closed)
        @test CSG._in(CartesianPoint(0.3, 0.0, -1.0), p_tube_closed)
        @test !CSG._in(CartesianPoint(0.0, 0.0, -1.5), p_tube_closed)
        @test !CSG._in(CartesianPoint(0.6, 0.6, 0.0), p_tube_closed)
        
        # PartialTopClosedTube (OpenPrimitive)
        p_tube_open   = CSG.PartialTopClosedTube{Float64, CSG.OpenPrimitive}(r_cone, phi_partial, hZ, origin, rotation)
        @test CSG._in(CartesianPoint(0.4, 0.1, -0.9), p_tube_open)
        @test !CSG._in(CartesianPoint(0.5, 0.0, -1.0), p_tube_open)
        @test !CSG._in(CartesianPoint(1.0, 0.0, 1.0), p_tube_open)
        @test !CSG._in(CartesianPoint(0.0, 0.0, -1.5), p_tube_open)
        @test !CSG._in(CartesianPoint(0.5, 0.0, -1.0), p_tube_open)
        @test !CSG._in(CartesianPoint(0.3, 0.0, -1.0), p_tube_open)
        
        surfs_closed = CSG.surfaces(p_tube_closed)
        @test surfs_closed[1] isa CSG.EllipticalSurface
        @test surfs_closed[2] isa CSG.PartialConeMantle
        @test surfs_closed[3] isa CSG.PartialConeMantle
        @test surfs_closed[4] isa CSG.Triangle
        @test surfs_closed[5] isa CSG.Triangle

        surfs_open = CSG.surfaces(p_tube_open)
        @test length(surfs_open) == 5
        @test surfs_open[1] isa CSG.EllipticalSurface
        @test surfs_open[2] isa CSG.PartialConeMantle
        @test surfs_open[3] isa CSG.PartialConeMantle
        @test surfs_open[4] isa CSG.Triangle
        @test surfs_open[5] isa CSG.Triangle

        # BottomClosedTube and PartialBottomClosedTube
        r_inner = 0.5
        r_outer = (0.8, 1.0)
        r_cone = (r_inner, r_outer)
        
        tube_closed       = CSG.BottomClosedTube{Float64, CSG.ClosedPrimitive}(r_cone, nothing, hZ, origin, rotation)
        tube_open         = CSG.BottomClosedTube{Float64, CSG.OpenPrimitive}(r_cone, nothing, hZ, origin, rotation)
        p_tube_closed     = CSG.PartialBottomClosedTube{Float64, CSG.ClosedPrimitive}(r_cone, phi_partial, hZ, origin, rotation)
        p_tube_open       = CSG.PartialBottomClosedTube{Float64, CSG.OpenPrimitive}(r_cone, phi_partial, hZ, origin, rotation)

        @test CSG._in(CartesianPoint(0.5, 0.0, -1.0), tube_closed)
        @test CSG._in(CartesianPoint(0.8, 0.0, 1.0), tube_closed)
        @test CSG._in(CartesianPoint(-0.9, 0.0, 0.9), tube_closed)
        @test !CSG._in(CartesianPoint(1.1, 0.0, 0.0), tube_closed)
        @test !CSG._in(CartesianPoint(0.5, 0.0, 1.1), tube_closed)
        
        @test !CSG._in(CartesianPoint(0.5, 0.0, -1.0), tube_open)
        @test !CSG._in(CartesianPoint(0.8, 0.0, 1.0), tube_open)
        @test CSG._in(CartesianPoint(-0.9, 0.0, 0.9), tube_open)
        @test !CSG._in(CartesianPoint(1.1, 0.0, 0.0), tube_open)
        @test !CSG._in(CartesianPoint(0.5, 0.0, 1.1), tube_open)

        @test CSG._in(CartesianPoint(0.5, 0.0, -1.0), p_tube_closed)
        @test CSG._in(CartesianPoint(0.8, 0.0, 1.0), p_tube_closed)
        @test CSG._in(CartesianPoint(1.0, 0.0, 1.0), p_tube_closed)
        @test CSG._in(CartesianPoint(0.5, 0.0, -1.0), p_tube_closed)
        @test !CSG._in(CartesianPoint(-0.9, 0.0, 0.9), p_tube_closed)
        @test CSG._in(CartesianPoint(0.9, 0.1, 0.9), p_tube_closed)
        @test !CSG._in(CartesianPoint(1.1, 0.0, 0.0), p_tube_closed)
        @test !CSG._in(CartesianPoint(0.5, 0.0, 1.1), p_tube_closed)

        @test !CSG._in(CartesianPoint(0.5, 0.0, -1.0), p_tube_open)
        @test !CSG._in(CartesianPoint(0.8, 0.0, 1.0), p_tube_open)
        @test !CSG._in(CartesianPoint(1.0, 0.0, 1.0), p_tube_open)
        @test !CSG._in(CartesianPoint(0.0, 0.0, -1.0), p_tube_open)
        @test !CSG._in(CartesianPoint(-0.9, 0.0, 0.9), p_tube_open)
        @test CSG._in(CartesianPoint(0.9, 0.1, 0.9), p_tube_open)
        @test !CSG._in(CartesianPoint(1.1, 0.0, 0.0), p_tube_open)
        @test !CSG._in(CartesianPoint(0.5, 0.0, 1.1), p_tube_open)
        
        surfs_closed = CSG.surfaces(tube_closed)
        @test length(surfs_closed) == 3
        @test surfs_closed[1] isa CSG.EllipticalSurface
        @test surfs_closed[2] isa CSG.FullConeMantle
        @test surfs_closed[3] isa CSG.FullConeMantle
        
        surfs_open = CSG.surfaces(tube_open)
        @test length(surfs_open) == 3
        @test surfs_open[1] isa CSG.EllipticalSurface
        @test surfs_open[2] isa CSG.FullConeMantle
        @test surfs_open[3] isa CSG.FullConeMantle
        
        surfs_p_closed = CSG.surfaces(p_tube_closed)
        @test length(surfs_p_closed) == 5
        @test surfs_p_closed[1] isa CSG.EllipticalSurface
        @test surfs_p_closed[2] isa CSG.PartialConeMantle
        @test surfs_p_closed[3] isa CSG.PartialConeMantle
        @test surfs_p_closed[4] isa CSG.Triangle
        @test surfs_p_closed[5] isa CSG.Triangle
        
        surfs_p_open = CSG.surfaces(p_tube_open)
        @test length(surfs_p_open) == 5
        @test surfs_p_open[1] isa CSG.EllipticalSurface
        @test surfs_p_open[2] isa CSG.PartialConeMantle
        @test surfs_p_open[3] isa CSG.PartialConeMantle
        @test surfs_p_open[4] isa CSG.Triangle
        @test surfs_p_open[5] isa CSG.Triangle
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
    
        # Test different constructor versions and dictionary construction
        torus1 = @inferred CSG.Torus(CSG.ClosedPrimitive, r_torus = 1f0, r_tube = 1f0)
        torus2 = @inferred CSG.Torus{Float32}(r_torus = 1.0, r_tube = 1.0)
        @test torus1 == torus2
        dict = Dict("torus"   => Dict(
                "r_torus"       => 1.0,
                "r_tube"        => 1.0))
                # "phi"           => π,
                # "theta"         => (π,2π)))
        torus4 = Geometry(T,dict,default_units,no_translations)
        output = Dictionary(torus4)
        @test dict == output

        ## Test in method for several tori
        # Non transformed geometries
        @test !in(CartesianPoint{Float32}(0,0,0.5),torus1)
        @test in(CartesianPoint{Float32}(0,0,0),torus1)
        torus1_open = @inferred CSG.Torus(CSG.OpenPrimitive, r_torus = 1f0, r_tube = 1f0)
        @test in(CartesianPoint{Float32}(0,2,0),torus1)
        @test !in(CartesianPoint{Float32}(0,2,0),torus1_open)
        # Boundary behaviour of closed and open primitives for transformed 
        # geometries (origin, rotation)
        torus_closed_trafo = @inferred CSG.Torus(CSG.ClosedPrimitive,r_torus=2.0, r_tube=3.0, φ = π, origin=CartesianPoint{Float32}(0,0,1),rotation=SMatrix{3}(0,0,-1,0,1,0,1,0,0))
        torus_open_trafo = @inferred CSG.Torus(CSG.OpenPrimitive,r_torus=2.0, r_tube=3.0, φ = π, origin=CartesianPoint{Float32}(0,0,1),rotation=SMatrix{3}(0,0,-1,0,1,0,1,0,0))

        @test in(CartesianPoint{Float64}(1,0,-1),torus_closed_trafo)
        @test !in(CartesianPoint{Float64}(1,0,-1),torus_open_trafo)
        @test in(CartesianPoint{Float64}(0,0,0),torus_closed_trafo)
        @test !in(CartesianPoint{Float64}(0,0,0),torus_open_trafo)
        tol = 1e-8
        @test !in(CartesianPoint{Float64}(0,-tol,0),torus_closed_trafo)
        @test in(CartesianPoint{Float64}(0,+tol,0),torus_open_trafo)
        #Test all valid keyword argument types
        @inferred CSG.Torus{Float64}(φ=nothing,θ=nothing)
        CSG.Torus{Float64}(φ=π,θ=(π,2π)) #inferred does not work here, due to evaluation of _get_conemantle_type(θ) during runtime
        #Test where phi is a Tuple
        tuple_torus = @inferred CSG.Torus{Float64}(φ=(π/4,3*π/4),θ=nothing,rotation=SMatrix{3}(0,0,-1,0,1,0,1,0,0))
        rot_torus = @inferred CSG.Torus{Float64}(φ=π/2,θ=nothing,rotation=SMatrix{3}(0,0,-1,0,1,0,1,0,0) * SMatrix{3}(cos(π/4),sin(π/4),0,-sin(π/4),cos(π/4),0,0,0,1))
        @test tuple_torus ==rot_torus

        @testset "_in" begin
        for r_tube in (2.0, (1.0, 2.0))
            r_torus = 5.0
            origin = CartesianPoint(0.0, 0.0, 0.0)
            rot = RotX(0.0)
            
            # Hollow OpenPrimitive Torus
            t_open = CSG.Torus{Float64, CSG.OpenPrimitive}(r_torus, r_tube, nothing, nothing, origin, rot)
            
            # Point safely inside tube (avoid boundaries due to csgtol)
            rmin, rmax = CSG._radial_endpoints(r_tube)
            safe_r = (rmin + rmax) / 2
            pt_inside = CartesianPoint(r_torus + safe_r, 0.0, 0.0)
            @test CSG._in(pt_inside, t_open)
            
            # Point outside radially   
            pt_out_rad = CartesianPoint(r_torus + rmax + 1.0, 0.0, 0.0)
            @test !CSG._in(pt_out_rad, t_open)

            # Hollow OpenPrimitive with φ restriction
            t_phi = CSG.Torus{Float64, CSG.OpenPrimitive}(r_torus, r_tube, (0.0, π/2), nothing, origin, rot)
            pt_out_phi = CartesianPoint(-r_torus - safe_r, 0.0, 0.0) # φ = π, outside
            @test !CSG._in(pt_out_phi, t_phi)

            t_phi = CSG.Torus{Float64, CSG.OpenPrimitive}(r_torus, r_tube, (0.0, 3π/2), nothing, origin, rot)
            pt_out_phi = CartesianPoint(-r_torus - safe_r, 0.0, 0.0) # φ = π, outside
            @test CSG._in(pt_out_phi, t_phi)

            # Hollow OpenPrimitive with θ restriction
            t_theta = CSG.Torus{Float64, CSG.OpenPrimitive}(r_torus, r_tube, nothing, (0.0, π/4), origin, rot)
            # Point with z above allowed θ
            pt_out_theta = CartesianPoint(r_torus + safe_r, 0.0, rmax + 1.0)
            @test !CSG._in(pt_out_theta, t_theta)

            t_theta = CSG.Torus{Float64, CSG.OpenPrimitive}(r_torus, r_tube, nothing, (0.0, 3π/2), origin, rot)
            # Point with z above allowed θ
            pt_out_theta = CartesianPoint(r_torus + safe_r, 0.0, rmax + 1.0)
            @test !CSG._in(pt_out_theta, t_theta)
        end
        end
        @testset "surfaces" begin
            r_torus = 5.0
            φ, θ = π/2, π/4
            origin = CartesianPoint(0.0, 0.0, 0.0)
            rot = RotX(0.0)

            # FullThetaTorus ClosedPrimitive
            r_tube = 1.0
            
            t_closed = CSG.Torus{Float64, CSG.ClosedPrimitive}(r_torus, r_tube, nothing, nothing, origin, rot)
            s = CSG.surfaces(t_closed)
            @test length(s) == 1

            t = CSG.Torus{Float64, CSG.ClosedPrimitive}(r_torus, r_tube, (0.0, φ), (0.0, θ), origin, rot)
            tm, es1, es2 = CSG.surfaces(t)

            @test tm isa CSG.TorusMantle{Float64,Float64,Tuple{Float64,Float64},:inwards}
            @test es1 isa CSG.EllipticalSurface
            @test es2 isa CSG.ConeMantle

            t_full_closed = CSG.FullThetaTorus{Float64, CSG.ClosedPrimitive}(r_torus, r_tube, φ, nothing, origin, rot)
            tm, es1, es2 = CSG.surfaces(t_full_closed)
            @test tm isa CSG.FullThetaTorusMantle
            @test es1 isa CSG.CircularArea
            @test es2 isa CSG.CircularArea
            
            # FullThetaTorus OpenPrimitive                                                 
            t_open = CSG.Torus{Float64, CSG.OpenPrimitive}(r_torus, r_tube, nothing, nothing, origin, rot)
            s_open = CSG.surfaces(t_open)
            @test length(s_open) == 1

            t = CSG.Torus{Float64, CSG.OpenPrimitive}(r_torus, r_tube, (0.0, φ), (0.0, θ), origin, rot)
            tm, es1, es2 = CSG.surfaces(t)
            
            @test tm isa CSG.TorusMantle{Float64,Float64,Tuple{Float64,Float64},:outwards}
            @test es1 isa CSG.EllipticalSurface
            @test es2 isa CSG.ConeMantle

            t_full_open   = CSG.FullThetaTorus{Float64, CSG.OpenPrimitive}(r_torus, r_tube, φ, nothing, origin, rot)
            tm, es1, es2 = CSG.surfaces(t_full_open)
            @test tm isa CSG.FullThetaTorusMantle
            @test es1 isa CSG.CircularArea
            @test es2 isa CSG.CircularArea
            
            # HollowThetaTorus ClosedPrimitive
            r_tube_hollow = (0.5, 1.0)
            
            t_hollow_closed = CSG.Torus{Float64, CSG.ClosedPrimitive}(r_torus, r_tube_hollow, (0.0, φ), (0.0, θ), origin, rot)
            s_hollow_closed = CSG.surfaces(t_hollow_closed)
            @test length(s_hollow_closed) == 6
            tm_in, tm_out, es1, es2 = CSG.surfaces(t_hollow_closed)
            
            @test tm_in  isa CSG.TorusMantle{Float64, Float64, Tuple{Float64, Float64}, :outwards}
            @test tm_out isa CSG.TorusMantle{Float64, Float64, Tuple{Float64, Float64}, :inwards}
            @test es1 isa CSG.EllipticalSurface
            @test es2 isa CSG.ConeMantle

            t_hollow_closed = CSG.HollowThetaTorus{Float64, CSG.ClosedPrimitive}(r_torus, r_tube_hollow, φ, nothing, origin, rot)
            tm_in, tm_out, es1, es2 = CSG.surfaces(t_hollow_closed)
            @test tm_in  isa CSG.FullThetaTorusMantle
            @test tm_out isa CSG.FullThetaTorusMantle
            @test es1    isa CSG.Annulus
            @test es2    isa CSG.Annulus
            
            # HollowThetaTorus OpenPrimitive
            t_hollow_open = CSG.Torus{Float64, CSG.OpenPrimitive}(r_torus, r_tube_hollow, (0.0, φ), (0.0, θ), origin, rot)
            s_hollow_open = CSG.surfaces(t_hollow_open)
            @test length(s_hollow_open) == 6
            tm_in, tm_out, es1, es2 = CSG.surfaces(t_hollow_open)

            @test tm_in  isa CSG.TorusMantle{Float64, Float64, Tuple{Float64, Float64}, :inwards}
            @test tm_out isa CSG.TorusMantle{Float64, Float64, Tuple{Float64, Float64}, :outwards}
            @test es1 isa CSG.EllipticalSurface
            @test es2 isa CSG.ConeMantle

            t_hollow_open = CSG.HollowThetaTorus{Float64, CSG.OpenPrimitive}(r_torus, r_tube_hollow, φ, nothing, origin, rot)
            tm_in, tm_out, es1, es2 = CSG.surfaces(t_hollow_open)
            @test tm_in  isa CSG.FullThetaTorusMantle
            @test tm_out isa CSG.FullThetaTorusMantle
            @test es1    isa CSG.Annulus
            @test es2    isa CSG.Annulus
        end
    end    
    @testset "TorusMantle" begin
        r_torus, r_tube = 5.0, 1.0
        origin = CartesianPoint{Float64}(0.0, 0.0, 0.0)
        rot = SMatrix{3,3,Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
        tm_out = CSG.TorusMantle{Float64}(:outwards; r_torus=r_torus, r_tube=r_tube, origin=origin, rotation=rot)
        tm_in  = CSG.TorusMantle{Float64}(:inwards; r_torus=r_torus, r_tube=r_tube, origin=origin, rotation=rot)
        
        # Test normal
        pt = CartesianPoint{Float64}(6.0, 0.0, 0.0)
        n_out = CSG.normal(tm_out, pt)
        n_in  = CSG.normal(tm_in, pt)
        
        @test n_out isa CartesianVector{Float64}
        @test n_in isa CartesianVector{Float64}
        @test dot(n_out, n_in) < -0.999
        
        # Test vertices
        verts = CSG.vertices(tm_out, 8)
        @test verts isa Vector{CartesianPoint{Float64}}
        @test length(verts) > 0
        
        # Test sample
        samples = CSG.sample(tm_out, 0.1)
        @test samples isa Vector{CartesianPoint{Float64}}
        @test length(samples) > 0
        
        # Test TorusThetaSurface
        r = 1.0
        φ = π/4
        hZ = 0.5
        
        ts_flat = CSG.TorusThetaSurface(r, φ, hZ, origin, rot, Val(:flat))
        @test ts_flat isa CSG.EllipticalSurface
        
        ts_in = CSG.TorusThetaSurface(r, φ, hZ, origin, rot, Val(:inwards))
        @test ts_in isa CSG.ConeMantle
        
        ts_out = CSG.TorusThetaSurface(r, φ, hZ, origin, rot, Val(:outwards))
        @test ts_out isa CSG.ConeMantle
    end
    @testset "Polycone" begin
        polycone1 = CSG.Polycone(CSG.ClosedPrimitive,r=Float32[0,2,4,0,0],z=Float32[0,1,2,3,0],origin=zero(CartesianPoint{Float16}),rotation=one(SMatrix{3, 3, Float16, 9}))
        polycone2 = CSG.Polycone{Float32}(r = [0,2,4,0], z = [0,1,2,3]) # omit last entry
        @test polycone1 == polycone2
        
        dict = Dict("polycone" => Dict(
            "r" => (0,2,4,0,0),
            "z" => (0,1,2,3,0),
        ))
        polycone3 = Geometry(T,dict,default_units,no_translations)
        @test dict == Dictionary(polycone3)

        polycone_open = CSG.Polycone(CSG.OpenPrimitive,r=Float32[0,2,4,0,0],z=Float32[0,1,2,3,0],origin=zero(CartesianPoint{Float16}),rotation=one(SMatrix{3, 3, Float16, 9}))
        @test  in(CartesianPoint{Float32}(4,0,2), polycone1)
        @test !in(CartesianPoint{Float32}(4,0,2), polycone_open)

        # Points should be sorted counter-clockwise (internally)
        dict2 = Dict("polycone" => Dict(
            "r" => (0,0,4,2,0),
            "z" => (0,3,2,1,0),
        ))
        polycone4 = Geometry(T,dict2,default_units,no_translations)
        @test polycone3 == polycone4

        # Extremum
        e = CSG.extremum(polycone3)
        @test all(e .>= abs.(polycone3.r)) && all(e .>= abs.(polycone3.z))
    end
    @testset "Ellipsoid" begin
        ellip1 = @inferred CSG.Ellipsoid(CSG.ClosedPrimitive,r=1f0,origin = zero(CartesianPoint{Float16}),rotation = one(SMatrix{3, 3, Float16, 9}))
        ellip2 = @inferred CSG.Ellipsoid{Float32}(r=1.0)
        @test ellip1 === ellip2
    
        dict = Dict("sphere" => Dict("r" => 1.0))
        ellip3 = Geometry(T,dict,default_units,no_translations)
        output = Dictionary(ellip3)
        @test dict == output

        # Partial ellipsoid not supported yet
        dict_φ = Dict("sphere" => Dict("r" => 1.0, "phi" => 3.0))
        @test_throws Exception Geometry(T,dict_φ,default_units,no_translations)
        dict_θ = Dict("sphere" => Dict("r" => 1.0, "theta" => 3.0))
        @test_throws Exception Geometry(T,dict_θ,default_units,no_translations)
    
        ellip_open = @inferred CSG.Ellipsoid(CSG.OpenPrimitive,r=1f0,origin = zero(CartesianPoint{Float16}),rotation = one(SMatrix{3, 3, Float16, 9}))
        @test in(CartesianPoint{Float32}(1,0,0),ellip1)
        @test !in(CartesianPoint{Float32}(1,0,0),ellip_open)
    
        ellip_closed_trafo = @inferred CSG.Ellipsoid(CSG.ClosedPrimitive,r=1.0, origin=CartesianPoint{Float64}(1/sqrt(3),1/sqrt(3),1/sqrt(3)),rotation=SMatrix{3}(0.5,sqrt(3)/2,0,-sqrt(3)/2,0.5,0,0,0,1))
        ellip_open_trafo = @inferred CSG.Ellipsoid(CSG.OpenPrimitive,r=1.0, origin=CartesianPoint{Float64}(1/sqrt(3),1/sqrt(3),1/sqrt(3)),rotation=SMatrix{3}(0.5,sqrt(3)/2,0,-sqrt(3)/2,0.5,0,0,0,1))
        @test !in(CartesianPoint{Float64}(-1e-8,0,0),ellip_closed_trafo)
        @test in(CartesianPoint{Float64}(1e-8,0,0),ellip_open_trafo)

        dict = Dict("sphere" => Dict("r" => Dict("from" => 1.0, "to" => 2.0)))
        ellip_hollow = Geometry(T,dict,default_units,no_translations)
        @test !in(CartesianPoint{T}(0,0,0),ellip_hollow)
        @test in(CartesianPoint{T}(1.0,0,0),ellip_hollow)
        @test in(CartesianPoint{T}(1.5,0,0),ellip_hollow)
        @test in(CartesianPoint{T}(2.0,0,0),ellip_hollow)

        # Extremum
        @test CSG.extremum(ellip1) >= ellip1.r
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

        # Test Volume Primitives
        pt_outside = CartesianPoint{Float32}(3, 0, 0)
        pt_inside  = CartesianPoint{Float32}(0.5, 0, 0)
        # Expected distances
        expected_distance_outside = 2.0f0  # nearest face along x
        expected_distance_inside  = 0.5f0  # nearest face along x
        @test CSG.distance(pt_outside, box2) ≈ expected_distance_outside
        @test CSG.distance(pt_inside, box2) ≈ expected_distance_inside
        # Extremum
        e = CSG.extremum(box1)
        @test e >= box1.hX && e >= box1.hY && e >= box1.hZ
    end
    @testset "RegularPrism" begin
        prism1 = @inferred CSG.RegularPrism{3}(CSG.ClosedPrimitive,r=1f0, hZ = 2f0)
        prism2 = @inferred CSG.RegularPrism{3,Float32}(r=1.0, hZ = 2f0)
        @test prism1 === prism2
        
        dict = Dict("TriangularPrism" => Dict("r" => 1.0, "h" => 1.0))
        prism = Geometry(T,dict,default_units,no_translations)
        output = Dictionary(prism)
        @test dict == output
        
        dict = Dict("QuadranglePrism" => Dict("r" => 1.0, "h" => 1.0))
        prism = Geometry(T,dict,default_units,no_translations)
        output = Dictionary(prism)
        @test dict == output

        dict = Dict("PentagonalPrism" => Dict("r" => 1.0, "h" => 1.0))
        prism = Geometry(T,dict,default_units,no_translations)
        output = Dictionary(prism)
        @test dict == output
        
        dict = Dict("HexagonalPrism"  => Dict("r" => 1.0, "h" => 1.0))
        prism = Geometry(T,dict,default_units,no_translations)
        output = Dictionary(prism)
        @test dict == output
        
        # Test in() method
        @test in(CartesianPoint{Float32}(0,0,0),prism1)
        @test !in(CartesianPoint{Float32}(1.5,0,0),prism1)
        tol = 1e-8
        prism_closed = @inferred CSG.RegularPrism{4}(CSG.ClosedPrimitive; r =1.,hZ =1., origin=CartesianPoint{Float32}(0,0,5),rotation=SMatrix{3}(0,0,-1,0,1,0,1,0,0))
        prism_open = @inferred CSG.RegularPrism{5}(CSG.OpenPrimitive; r =1.,hZ =1., origin=CartesianPoint{Float32}(0,0,5),rotation=SMatrix{3}(0,0,-1,0,1,0,1,0,0))
        @test in(CartesianPoint{Float64}(0,0,4),prism_closed)
        @test !in(CartesianPoint{Float64}(0,0,4),prism_open)
        @test !in(CartesianPoint{Float64}(0,0,4-tol),prism_closed)
        @test in(CartesianPoint{Float64}(0,0,4+tol),prism_open)
        
        #Test different prism
        @inferred CSG.RegularPrism{3,Float32}()
        @inferred CSG.RegularPrism{4,Float32}()
        @inferred CSG.RegularPrism{5,Float32}()
        @inferred CSG.RegularPrism{6,Float32}()

        # Extremum
        e = CSG.extremum(prism_closed)
        @test e >= prism_closed.r && e >= prism_closed.hZ
    end
    @testset "EllipsoidMantle" begin
        em = CSG.EllipsoidMantle{Float32}(r = 1f0)
        CSG.EllipsoidMantle() 

        # Point on the surface
        pt = CartesianPoint{Float32}(1f0, 0f0, 0f0)

        # Test normal
        normal_vec = CSG.normal(em, pt)
        @test normal_vec isa CartesianVector{Float32}
        @test abs(norm(normal_vec) - 1f0) < 1e-5  # normalized
        
        em_out = CSG.EllipsoidMantle{Float32}(:outwards; r=(1f0, 2f0, 3f0))
        em_in  = CSG.EllipsoidMantle{Float32}(:inwards;  r=(1f0, 2f0, 3f0))
        pt2 = CartesianPoint{Float32}(0.6f0, 0.4f0, 0.2f0)

        n_out = CSG.normal(em_out, pt2)
        n_in  = CSG.normal(em_in,  pt2)

        @test n_out isa CartesianVector{Float32}
        @test n_in isa CartesianVector{Float32}
        v_pt2 = CartesianVector(pt2.x, pt2.y, pt2.z)
        # Outward normal must point away from center
        @test dot(n_out, v_pt2) > 0
        # Inward normal must point toward center
        @test dot(n_in,  v_pt2) < 0
        
        p_local = CSG._transform_into_object_coordinate_system(pt2, em_out)
        x = p_local.x
        y = p_local.y
        z = p_local.z

        expected_local = CartesianVector(
            sign(x) * (x / em_out.r[1])^2,
            sign(y) * (y / em_out.r[2])^2,
            sign(z) * (z / em_out.r[3])^2
        )

        expected_global = CSG._transform_into_global_coordinate_system(expected_local, em_out)

        @test norm(n_out - expected_global) < 1e-6
        @test norm(n_in  + expected_global) < 1e-6
        
        # Test vertices
        verts = CSG.vertices(em, 8)
        @test verts isa Vector{CartesianPoint{Float32}}
        @test length(verts) > 0
	
	# Test sample
        samples = CSG.sample(em, 0.1f0)
        @test samples isa Vector{CartesianPoint{Float32}}
        @test length(samples) > 0
    end
    @testset "EllipticalSurface" begin
        @inferred CSG.EllipticalSurface{Float32}(r = 1f0)
        @inferred CSG.EllipticalSurface()
        @inferred CSG.EllipticalSurface(φ=(1.,2.))    

        es_full    = CSG.EllipticalSurface{Float32}(r=1f0)
        es_partial = CSG.EllipticalSurface{Float32}(r=1f0, φ=Float32(π/2))
        es_annulus = CSG.EllipticalSurface{Float32}(r=(0.5f0, 1f0))

        verts_full = CSG.vertices(es_full, 8)
        @test verts_full isa Vector{CartesianPoint{Float32}}
        @test length(verts_full) == 10

        verts_annulus = CSG.vertices(es_annulus, 8)
        @test verts_annulus isa Vector{CartesianPoint{Float32}}
        @test length(verts_annulus) == 2*9

        samples_full = CSG.sample(es_full, 0.1f0)
        @test samples_full isa Vector{CartesianPoint{Float32}}
        @test length(samples_full) > 0

        samples_annulus = CSG.sample(es_annulus, 0.1f0)
        @test samples_annulus isa Vector{CartesianPoint{Float32}}
        @test length(samples_annulus) > 0

        circ, edges = CSG.lines(CSG.EllipticalSurface(r=1f0, φ=Float32(π/2)))
        @test circ isa CSG.Ellipse{Float32, Float32, Float32}
        @test all(isa.(edges, CSG.Edge{Float32}))

        pt = CartesianPoint{Float32}(1f0, 0f0, 0f0)
        dist_full = CSG.distance_to_surface(pt, es_full)
        @test dist_full isa Float32
        @test dist_full ≈ 0f0 atol=1e-5

        pt_off = CartesianPoint{Float32}(1.5f0, 0f0, 0f0)
        dist_annulus = CSG.distance_to_surface(pt_off, es_annulus)
        @test dist_annulus isa Float32
        @test dist_annulus > 0f0
    end

    @testset "distance_to_surface" begin
        # Test distance_to_surface and _above_or_below_polygon accept AbstractCoordinatePoint
        cm = CSG.ConeMantle{T}()
        es = CSG.EllipticalSurface{T}()
        torus = CSG.TorusMantle{T}()
        pts = SVector(
            CartesianPoint{T}(0.0, 0.0, 0.0),
            CartesianPoint{T}(1.0, 0.0, 0.0),
            CartesianPoint{T}(1.0, 1.0, 0.0),
            CartesianPoint{T}(0.0, 1.0, 0.0),
            CartesianPoint{T}(0.5, 0.1, 0.0)
        )
        quad = CSG.Polygon{4,T}(pts[1:4])
        poly = CSG.Polygon{5,T}(pts)
        pt = CartesianPoint{T}(0.1, 0.1, 0.1)

        @test CSG.distance_to_surface(pt, cm) isa T
        @test CSG.distance_to_surface(pt, es) isa T
        @test CSG.distance_to_surface(pt, torus) isa T
        @test CSG.distance_to_surface(pt, quad) isa T
        @test CSG.distance_to_surface(pt, poly) isa T
        @test_throws ArgumentError CSG.distance_to_surface(pt, CSG.CSGDifference(CSG.Box{T}(), CSG.Cone{T}()))
    end
    @testset "Plane" begin
        plane1 = @inferred CSG.Plane{Float32}(normal = CartesianVector(1,0,0))
        plane2 = @inferred CSG.Plane(normal = CartesianVector(1f0,0f0,0f0))
        @test plane1 == plane2 

        # Test Plane functions use vector differences correctly
        p = CSG.Plane(CartesianPoint(0.0, 0.0, 0.0), CartesianVector(0.0, 0.0, 1.0))
        pt_front = CartesianPoint(0.0, 0.0, 1.0)
        diff = pt_front - CSG.origin(p)
        @test diff isa CartesianVector
        @test @inferred(CSG.isinfront(pt_front, p)) === true
        @test @inferred(CSG.isbehind(pt_front, p)) === false
        @test @inferred(CSG._distance(pt_front, p)) == 1.0
        @test @inferred(CSG.distance(pt_front, p)) == 1.0
    end
    @testset "ConeMantle" begin
        CSG.ConeMantle{Float32}(r = 1f0)
        CSG.ConeMantle()
        CSG.ConeMantle(φ=(1.,2.))    

        cm_out = CSG.ConeMantle{Float32}( :outwards; r=(1f0, 2f0), φ=nothing, hZ=1f0 )
        cm_in  = CSG.ConeMantle{Float32}( :inwards;  r=(1f0, 2f0), φ=nothing, hZ=1f0 )
        pt     = CartesianPoint{Float32}(0.7f0, 0.4f0, 0.0f0)

        normal_out = CSG.normal(cm_out, pt)
        normal_in  = CSG.normal(cm_in, pt)

        @test normal_out isa CartesianVector{Float32}
        @test normal_in isa CartesianVector{Float32}
        # inwards vs outwards should point opposite
        @test dot(normal_out, normal_in) ≈ -1f0 atol=1e-5

        v_pt = CartesianVector(pt.x, pt.y, pt.z)
        @test dot(normal_out, v_pt) > 0
        @test dot(normal_in,  v_pt) < 0

        # Cylindrical ConeMantle
        cm_cyl_out = CSG.ConeMantle{Float32}(:outwards; r=1f0, φ=nothing, hZ=1f0)
        cm_cyl_in = CSG.ConeMantle{Float32}(:inwards; r=1f0, φ=nothing, hZ=1f0)

        pt_cyl = CartesianPoint{Float32}(0.6f0, 0.8f0, 0.2f0)

        n_out = CSG.normal(cm_cyl_out, pt_cyl)
        n_in  = CSG.normal(cm_cyl_in,  pt_cyl)

        @test n_out isa CartesianVector{Float32}
        @test n_in  isa CartesianVector{Float32}

        # purely radial, no z component
        @test abs(n_out.z) ≤ 1f-6
        @test abs(n_in.z)  ≤ 1f-6
        @test dot(n_out, n_in) ≈ -1f0 atol=1f-6

        # radial direction test
        v_radial = normalize(CartesianVector(pt_cyl.x, pt_cyl.y, 0f0))
        @test dot(n_out, v_radial) ≈  1f0 atol=1f-6
        @test dot(n_in,  v_radial) ≈ -1f0 atol=1f-6
        
        n_arc = 8
        verts = CSG.vertices(cm_out, n_arc)
        @test verts isa Vector{CartesianPoint{Float32}}
        @test length(verts) == 2*(n_arc + 1)   # bottom + top circle points
        
        # All vertices should lie approximately on their respective radius circles
        r_bot = CSG.radius_at_z(cm_out, -cm_out.hZ)
        r_top = CSG.radius_at_z(cm_out, cm_out.hZ)
        for (i, v) in enumerate(verts)
            if i <= n_arc+1
                @test abs(norm(CartesianVector(v.x, v.y, 0)) - r_bot) < 1e-5
            else
                @test abs(norm(CartesianVector(v.x, v.y, 0)) - r_top) < 1e-5
            end
        end
        
        pts = CSG.sample(cm_out, Float32(1))
        @test pts isa Vector{CartesianPoint{Float32}}
        @test length(pts) >= 2  # must have at least 2 points

        # Check that all points are within height bounds
        for pt in pts
            @test -cm_out.hZ - 1e-5 ≤ pt.z ≤ cm_out.hZ + 1e-5
        end
    end
    @testset "Ellipse" begin
        ell1 = @inferred CSG.Ellipse{Float32}(r = 1.0f0, φ=10.0f0)
        ell2 = @inferred CSG.Ellipse(r = 1.0f0, φ=10.0f0) 
        @test ell1 == ell2
        @test CSG.extremum(ell1) >= ell1.r
        ell3 = @inferred CSG.Ellipse{T}(r = (1,2))
        @test all(CSG.extremum(ell3) .>= ell3.r)
    end
    @testset "Line" begin
        line1 = @inferred CSG.Line{Float32}()
        line2 = @inferred CSG.Line(direction = CartesianVector{Float32}(0,0,1)) 
        @test line1 == line2
        
        # Line along the x-axis, origin at (0,0,0)
        a = CartesianPoint{T}(0,0,0)
        b = CartesianPoint{T}(1,0,0)
        l = CSG.Line(a, b)
        
        # Point on the line → distance = 0
        pt_on_line = CartesianPoint{T}(2,0,0)
        @test CSG.distance(pt_on_line, l) ≈ T(0)

        # Point above the line in z direction → distance = 1
        pt_above = CartesianPoint{T}(2,0,1)
        @test CSG.distance(pt_above, l) ≈ T(1)
        
        # Point above the line in y direction → distance = 2
        pt_side = CartesianPoint{T}(3,2,0)
        @test CSG.distance(pt_side, l) ≈ T(2)
        
        # Line along arbitrary direction
        a2 = CartesianPoint{T}(1,1,1)
        b2 = CartesianPoint{T}(2,3,2)
        l2 = CSG.Line(a2, b2)
        pt2 = CartesianPoint{T}(3,1,1)
        # Compute expected distance manually
        v = pt2 - l2.origin
        d_expected = norm(cross(v, l2.direction)) / norm(l2.direction)
        @test CSG.distance(pt2, l2) ≈ d_expected 
    end
    @testset "Edge" begin
        edge1 = @inferred CSG.Edge{Float32}()
        edge2 = @inferred CSG.Edge(a = CartesianPoint{Float32}(0,0,0), b = CartesianPoint{Float16}(0,0,1)) 
        @test edge1 == edge2
        
        # Edge along x-axis
        a = CartesianPoint{T}(0,0,0)
        b = CartesianPoint{T}(1,0,0)
        e = CSG.Edge(a, b)
        
        # Case 1: Point on the segment → distance = 0
        pt_on = CartesianPoint{T}(0.5,0,0)
        @test CSG.distance(pt_on, e) ≈ 0.0
        pt_on = CartesianPoint{T}(0.3,0,0)
        @test CSG.distance_to_line(pt_on, e) ≈ 0.0
        
        # Case 2: Point beyond b → distance = ||pt - b||
        pt_beyond = CartesianPoint{T}(2,1,0)
        @test CSG.distance(pt_beyond, e) ≈ norm(pt_beyond - b)
        pt_beyond = CartesianPoint{T}(2,0,0)
        @test CSG.distance_to_line(pt_beyond, e) ≈ 1.0
        
        # Case 3: Point before a → distance = ||pt - a||
        pt_before = CartesianPoint{T}(-1,1,0)
        @test CSG.distance(pt_before, e) ≈ norm(pt_before - a)
        pt_before = CartesianPoint{T}(-1,0,0)
        @test CSG.distance_to_line(pt_before, e) ≈ 1.0
        
        # Case 4: Point off the segment but projection falls on segment
        pt_off = CartesianPoint{T}(0.5,2,0)
        v = b - a
        expected = norm(cross(pt_off - a, v)) / norm(v)
        @test CSG.distance(pt_off, e) ≈ expected
        pt_inside = CartesianPoint{T}(0.5,1,0)
        @test CSG.distance_to_line(pt_inside, e) ≈ 1.0

        e_vert = CSG.Edge(CartesianPoint{T}(0,0,0), CartesianPoint{T}(0,0,1))
        pt_off_vert = CartesianPoint{T}(1,0,0.5)
        @test CSG.distance(pt_off_vert, e_vert) ≈ T(1) 
    end
    @testset "Polygon" begin
        # Simple square in xy-plane
        pts = SVector(
            CartesianPoint{T}(0,0,0),
            CartesianPoint{T}(1,0,0),
            CartesianPoint{T}(1,1,0),
            CartesianPoint{T}(0,1,0)
        )
        poly = CSG.Polygon(pts)
        tri = CSG.Triangle{T}([pts[1], pts[2], pts[3]])
        
        @testset "_sample_excluding_border" begin
            samples = CSG._sample_excluding_border(tri, 0.2)
            @test all(p -> p isa CartesianPoint{T}, samples)
            @test all(p -> p[3] == 0, samples)  # xy-plane
            @test length(samples) > 0
        end
        
        @testset "sample" begin
            samples = CSG.sample(poly, 0.3)
            @test all(p -> p isa CartesianPoint{T}, samples)
            @test all(p -> p[3] == 0, samples)
            @test length(samples) > 0
        end

        @testset "distance" begin
            p_inside  = CartesianPoint{T}(0.5,0.5,0)
            p_above   = CartesianPoint{T}(0.5,0.5,1)
            p_outside = CartesianPoint{T}(2,0.5,0)
            
            @test isapprox(CSG.distance(p_inside, poly), 0; atol=1e-12)
            @test isapprox(CSG.distance(p_above,  poly), 1; atol=1e-12)
            @test isapprox(CSG.distance(p_outside, poly), 1; atol=1e-12)
        end

        @testset "_get_rot_for_rotation_on_xy_plane" begin
        poly_xy = CSG.Polygon(SVector(
            CartesianPoint{T}(0,0,0),
            CartesianPoint{T}(1,0,0),
            CartesianPoint{T}(1,1,0),
            CartesianPoint{T}(0,1,0),
            ))
            
        rot_xy = CSG._get_rot_for_rotation_on_xy_plane(poly_xy)
            
            pts_tilted = SVector(
                CartesianPoint{T}(0,0,0),
                CartesianPoint{T}(1,0,1),
                CartesianPoint{T}(1,1,1),
                CartesianPoint{T}(0,1,0),
            )
            poly_tilted = CSG.Polygon(pts_tilted)
            
            rot = CSG._get_rot_for_rotation_on_xy_plane(poly_tilted)
            # Rotated vertices should lie on z = 0
            for v in (rot * (p - cartesian_zero) for p in pts_tilted)
                @test abs(v[3]) ≤ 1e-12
            end
            
            # Rotated normal should be (0,0,1)
            n_rot = rot * SVector{3,T}(CSG.normal(poly_tilted)...)
            @test isapprox(n_rot, SVector(0,0,1), atol=1e-12)
        end
       
        @testset "in polygon" begin
            # Inside point
            p_in = CartesianPoint{T}(0.5, 0.5, 0)
            @test CSG.in(p_in, poly)
            
            # Outside point
            p_out = CartesianPoint{T}(1.5, 0.5, 0)
            @test !CSG.in(p_out, poly)
            
            # On boundary
            p_edge = CartesianPoint{T}(0.5, 0, 0)
            @test CSG.in(p_edge, poly)   # on-edge should return true
            
            # Off the plane
            p_above_plane = CartesianPoint{T}(0.5, 0.5, 1)
            @test !CSG.in(p_above_plane, poly)

            # Tilted Polygon
            pts_tilted = SVector(
                CartesianPoint{T}(0,0,0),
                CartesianPoint{T}(1,0,1),
                CartesianPoint{T}(1,1,1),
                CartesianPoint{T}(0,1,0)
            )
            poly_tilt = CSG.Polygon(pts_tilted)
            
            # Inside projected area
            p_in_tilted = CartesianPoint{T}(0.5, 0.5, 0.5)
            @test CSG.in(p_in_tilted, poly_tilt)
            
            # Outside projected area
            p_out_tilted = CartesianPoint{T}(2, 0.5, 1)
            @test !CSG.in(p_out_tilted, poly_tilt)
        end
    end
end

@testset "Test geometry parsing" begin
    sim = Simulation{T}(SSD_examples[:InvertedCoax])
    c = sim.detector.contacts[2].geometry
    @test Geometry(T, Dictionary(c), default_units, no_translations) == c
end

@testset "Test geometry rounding" begin
    pt = CartesianPoint{Float64}(1.23456789, 0.0, 1e-13)
    rounded_pt = CSG.geom_round(pt)
    
    @test rounded_pt.x == 1.23456789
    @test rounded_pt.y == 0.0
    @test rounded_pt.z == 0.0
    
    pt32 = CartesianPoint{Float32}(1.2345678f0, Float32(1e-7), Float32(-1e-7))
    rounded_pt32 = CSG.geom_round(pt32)
    
    @test rounded_pt32.x ≈ 1.234568f0
    @test rounded_pt32.y == 0.0f0
    @test rounded_pt32.z == 0.0f0
        
    pts = [CartesianPoint{Float64}(1e-13, 2.3456789, 0.0),
           CartesianPoint{Float64}(0.123456789, 1e-14, -0.987654321)]
    rounded_pts = CSG.geom_round(pts)
    
    @test rounded_pts[1].x == 0.0
    @test rounded_pts[1].y ≈ 2.3456789
    @test rounded_pts[1].z == 0.0
    @test rounded_pts[2].x ≈ 0.123456789
    @test rounded_pts[2].y == 0.0
    @test rounded_pts[2].z ≈ -0.987654321
end

@testset "Geometry Test" begin
    outer_transformations = (
        rotation    = SMatrix{3,3,T}(I),
        translation = CartesianVector{T}(1.0, 2.0, 3.0)
    )
    dict = Dict(
        "Z" => π/2,
        "box" => Dict(
            "hX"        => 1.0,
            "hY"        => 2.0,
            "hZ"        => 3.0,
            "translation" => CartesianVector{T}(4.0, 5.0, 6.0)
        )
    )
    input_units = (angle  = u"rad", length = u"m")
    geom = CSG.Geometry(T, CSG.Rotation, dict, input_units, outer_transformations)

    primitive_translation_rotated = geom.rotation * dict["box"]["translation"]

    expected_origin = CartesianPoint(
        outer_transformations.translation[1] + primitive_translation_rotated[1],
        outer_transformations.translation[2] + primitive_translation_rotated[2],
        outer_transformations.translation[3] + primitive_translation_rotated[3]
    )
    @test geom.origin == expected_origin

    expected_rot = SMatrix{3,3,T,9}(RotZ(π/2))
    @test geom.rotation ≈ expected_rot

    # --- Test + operator for AbstractConstructiveGeometry ---
    v = CartesianVector{T}(1.0, -2.0, 0.5)
    geom_translated = geom + v
    @test geom_translated.origin == geom.origin + v
    @test geom_translated.rotation ≈ geom.rotation
    # --- Test + operator for AbstractPrimitive ---
    prim = CSG.Box(CSG.ClosedPrimitive, hX=1.0, hY=2.0, hZ=3.0)
    prim_translated = prim + v
    @test prim_translated.origin == prim.origin + v
    @test prim_translated.rotation ≈ prim.rotation
    # --- Test * operator for AbstractPrimitive ---
    r = RotXYZ{T}(1.0, 1.0, 1.0)
    prim_rotated = r * prim
    @test prim_rotated.rotation ≈ r
    @test haskey(Dictionary(prim_rotated)["box"]["rotation"], "M")
end

@testset "_in_φ" begin
    # Cartesian points
    p1 = CartesianPoint(1.0, 0.0, 0.0)   # φ = 0
    p2 = CartesianPoint(0.0, 1.0, 0.0)   # φ = π/2
    p3 = CartesianPoint(-1.0, 0.0, 0.0)  # φ = π
    p4 = CartesianPoint(0.0, -1.0, 0.0)  # φ = 3π/2
    
    @test CSG._in_φ(p1, nothing)
    @test CSG._in_φ(p2, nothing)
    @test CSG._in_φ(p3, nothing)
    @test CSG._in_φ(p4, nothing)

    @test CSG._in_φ(p1, π/4)
    @test !CSG._in_φ(p2, π/4)
    @test CSG._in_φ(p2, π/2)
    @test !CSG._in_φ(p3, π/2)

    # Cylindrical points
    c1 = CylindricalPoint(1.0, 0.0, 0.0)      # φ = 0
    c2 = CylindricalPoint(1.0, π/2, 0.0)      # φ = π/2
    c3 = CylindricalPoint(1.0, π, 0.0)        # φ = π
    c4 = CylindricalPoint(1.0, 3π/2, 0.0)     # φ = 3π/2
    
    @test CSG._in_φ(c1, nothing)
    @test CSG._in_φ(c2, nothing)
    @test CSG._in_φ(c3, nothing)
    @test CSG._in_φ(c4, nothing)

    @test CSG._in_φ(c1, π/4)
    @test !CSG._in_φ(c2, π/4)
    @test CSG._in_φ(c2, π/2)
    @test !CSG._in_φ(c3, π/2)
end
    
@testset "CSG intersection" begin
    # Test in() accepts AbstractCoordinatePoint
    sphere = CSG.Ellipsoid(r = 1.0)
    box = CSG.Box(hX = 1.0, hY = 1.0, hZ = 1.0)
    inter = CSG.CSGIntersection(sphere, box)
    pt_in = CartesianPoint(0.2, 0.2, 0.2)
    @test pt_in in inter
end

@testset "distance_to_surface on segmented surfaces" begin
    pt_inφ1  = CartesianPoint{T}(0.0,  3.0, 0.0)
    pt_inφ2  = CartesianPoint{T}(0.0,  2.0, 0.0)
    pt_outφ1 = CartesianPoint{T}(2.0, -1.0, 0.0)
    pt_outφ2 = CartesianPoint{T}(0.0, -2.0, 0.0)


    # Segmented EllipticalSurface
    seg_ellipse = CSG.EllipticalSurface{T}(r = (1.0,2.0), origin = CartesianPoint(0.0, 0.0, 0.0), φ = (0.0, π/2))
    @test CSG.distance_to_surface(pt_inφ1,  seg_ellipse) ≈ T(1.0)
    @test CSG.distance_to_surface(pt_inφ2,  seg_ellipse) ≈ T(0.0)
    @test CSG.distance_to_surface(pt_outφ1, seg_ellipse) ≈ T(1.0)
    @test CSG.distance_to_surface(pt_outφ2, seg_ellipse) ≈ T(sqrt(5))

    # Segmented ConeMantle
    seg_cone = CSG.ConeMantle{T}(:outwards; r = (2.0, 2.0), φ = (0.0, π/2), hZ = 1.0)
    @test CSG.distance_to_surface(pt_inφ1,  seg_cone) ≈ T(1.0)
    @test CSG.distance_to_surface(pt_inφ2,  seg_cone) ≈ T(0.0)
    @test CSG.distance_to_surface(pt_outφ1, seg_cone) ≈ T(1.0)
    @test CSG.distance_to_surface(pt_outφ2, seg_cone) ≈ T(sqrt(8))

    @test CSG.distance_to_surface(pt_inφ1 + CartesianVector{T}(0,0,1.0),  seg_cone) ≈ T(1.0)
    @test CSG.distance_to_surface(pt_inφ1 - CartesianVector{T}(0,0,1.0),  seg_cone) ≈ T(1.0)
    @test CSG.distance_to_surface(pt_inφ2 + CartesianVector{T}(0,0,1.0),  seg_cone) ≈ T(0.0)
    @test CSG.distance_to_surface(pt_inφ2 - CartesianVector{T}(0,0,1.0),  seg_cone) ≈ T(0.0)
    @test CSG.distance_to_surface(pt_outφ1 + CartesianVector{T}(0,0,1.0), seg_cone) ≈ T(1.0)
    @test CSG.distance_to_surface(pt_outφ1 - CartesianVector{T}(0,0,1.0), seg_cone) ≈ T(1.0)
    @test CSG.distance_to_surface(pt_outφ2 + CartesianVector{T}(0,0,1.0), seg_cone) ≈ T(sqrt(8))
    @test CSG.distance_to_surface(pt_outφ2 - CartesianVector{T}(0,0,1.0), seg_cone) ≈ T(sqrt(8))
end
