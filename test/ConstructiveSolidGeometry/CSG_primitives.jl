using Test
using SolidStateDetectors
using StaticArrays
using LinearAlgebra
using Rotations

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
            r_torus, r_tube = 5.0, (1.0, 2.0)
            origin = CartesianPoint(0.0, 0.0, 0.0)
            rot = RotX(0.0)
            
            # Hollow OpenPrimitive Torus
            t_open = CSG.Torus{Float64, CSG.OpenPrimitive}(r_torus, r_tube, nothing, nothing, origin, rot)
            
            # Point safely inside tube (avoid boundaries due to csgtol)
            rmin, rmax = r_tube
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

            # Hollow OpenPrimitive with θ restriction
            t_theta = CSG.Torus{Float64, CSG.OpenPrimitive}(r_torus, r_tube, nothing, (0.0, π/4), origin, rot)
            # Point with z above allowed θ
            pt_out_theta = CartesianPoint(r_torus + safe_r, 0.0, rmax + 1.0)
            @test !CSG._in(pt_out_theta, t_theta)
        end
        @testset "surfaces" begin
            r_torus, r_tube = 5.0, 1.0
	    φ, θ = π/2, π/4
            origin = CartesianPoint(0.0, 0.0, 0.0)
            rot = RotX(0.0)

            # FullThetaTorus ClosedPrimitive
            t_closed = CSG.Torus{Float64, CSG.ClosedPrimitive}(r_torus, r_tube, nothing, nothing, origin, rot)
            s = CSG.surfaces(t_closed)
            @test length(s) == 1

            # FullThetaTorus OpenPrimitive                                                 
            t_open = CSG.Torus{Float64, CSG.OpenPrimitive}(r_torus, r_tube, nothing, nothing, origin, rot)
            s_open = CSG.surfaces(t_open)
            @test length(s_open) == 1
            # HollowThetaTorus ClosedPrimitive
            t_hollow = CSG.Torus{Float64, CSG.ClosedPrimitive}(r_torus, (0.5,1.0), (0.0, φ), (0.0, θ), origin, rot)
            s_hollow = CSG.surfaces(t_hollow)
            @test length(s_hollow) == 6
            
            # HollowThetaTorus OpenPrimitive
            t_hollow_open = CSG.Torus{Float64, CSG.OpenPrimitive}(r_torus, (0.5,1.0), (0.0, φ), (0.0, θ), origin, rot)
            s_hollow_open = CSG.surfaces(t_hollow_open)
            @test length(s_hollow_open) == 6
        end
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
    
        ellip_closed_trafo = @inferred CSG.Ellipsoid(CSG.ClosedPrimitive,r=1.0, origin=CartesianPoint{Float64}(1/sqrt(3),1/sqrt(3),1/sqrt(3)),rotation=SMatrix{3}(0.5,sqrt(3)/2,0,-sqrt(3)/2,0.5,0,0,0,1))
        ellip_open_trafo = @inferred CSG.Ellipsoid(CSG.OpenPrimitive,r=1.0, origin=CartesianPoint{Float64}(1/sqrt(3),1/sqrt(3),1/sqrt(3)),rotation=SMatrix{3}(0.5,sqrt(3)/2,0,-sqrt(3)/2,0.5,0,0,0,1))
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
    @testset "RegularPrism" begin
        prism1 = @inferred CSG.RegularPrism{3}(CSG.ClosedPrimitive,r=1f0, hZ = 2f0)
        prism2 = @inferred CSG.RegularPrism{3,Float32}(r=1.0, hZ = 2f0)
        @test prism1 === prism2
        
        dict = Dict("TriangularPrism"   => Dict(
                "r"       => 1.0,
                "h"        => 1.0))
        prism = Geometry(T,dict,default_units,no_translations)
        output = Dictionary(prism)
        @test dict == output
        
        dict = Dict("QuadranglePrism"   => Dict(
                "r"       => 1.0,
                "h"      =>1.0))
        prism = Geometry(T,dict,default_units,no_translations)
        output = Dictionary(prism)
        @test dict == output
        
        dict = Dict("HexagonalPrism"   => Dict(
                "r"       => 1.0,
                "h"      =>1.0))
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
    end
    @testset "EllipsoidMantle" begin
        em = CSG.EllipsoidMantle{Float32}(r = 1f0)
        CSG.EllipsoidMantle() 

        # Point on the surface
        pt = CartesianPoint{Float32}(1f0, 0f0, 0f0)

        # Test normal
        normal_vec = CSG.normal(em, pt)
        @test isa(normal_vec, CartesianVector{Float32})
        @test abs(norm(normal_vec) - 1f0) < 1e-5  # normalized
        
        # Test vertices
        verts = CSG.vertices(em, 8)
        @test isa(verts, Vector{CartesianPoint{Float32}})
        @test length(verts) > 0
	
	# Test sample
        samples = CSG.sample(em, 0.1f0)
        @test isa(samples, Vector{CartesianPoint{Float32}})
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
        @test isa(verts_full, Vector{CartesianPoint{Float32}})
        @test length(verts_full) == 10

        verts_annulus = CSG.vertices(es_annulus, 8)
        @test isa(verts_annulus, Vector{CartesianPoint{Float32}})
        @test length(verts_annulus) == 2*9

        samples_full = CSG.sample(es_full, 0.1f0)
        @test isa(samples_full, Vector{CartesianPoint{Float32}})
        @test length(samples_full) > 0

        samples_annulus = CSG.sample(es_annulus, 0.1f0)
        @test isa(samples_annulus, Vector{CartesianPoint{Float32}})
        @test length(samples_annulus) > 0

        circ, edges = CSG.lines(CSG.EllipticalSurface(r=1f0, φ=Float32(π/2)))
        @test isa(circ, CSG.Ellipse{Float32, Float32, Float32})
        @test all(isa.(edges, CSG.Edge{Float32}))

        p_inside = CartesianPoint{Float32}(0.5f0, 0f0, 0f0)
        p_outside = CartesianPoint{Float32}(1.5f0, 0f0, 0f0)
        @test CSG._in_cyl_r(p_inside, 1f0)
        @test !CSG._in_cyl_r(p_outside, 1f0)
        @test CSG._in_cyl_r(p_inside, (0.2f0, 0.6f0))
        @test !CSG._in_cyl_r(p_outside, (0.2f0, 0.6f0))

        pt = CartesianPoint{Float32}(1f0, 0f0, 0f0)
        dist_full = CSG.distance_to_surface(pt, es_full)
        @test isa(dist_full, Float32)
        @test dist_full ≈ 0f0 atol=1e-5

        pt_off = CartesianPoint{Float32}(1.5f0, 0f0, 0f0)
        dist_annulus = CSG.distance_to_surface(pt_off, es_annulus)
        @test isa(dist_annulus, Float32)
        @test dist_annulus > 0f0
    end
    @testset "Plane" begin
        plane1 = @inferred CSG.Plane{Float32}(normal = CartesianVector(1,0,0))
        plane2 = @inferred CSG.Plane(normal = CartesianVector(1f0,0f0,0f0))
        @test plane1 == plane2 
    end
    @testset "ConeMantle" begin
        CSG.ConeMantle{Float32}(r = 1f0)
        CSG.ConeMantle()
        CSG.ConeMantle(φ=(1.,2.))    

        cm_out = CSG.ConeMantle{Float32}( :outwards; r=(1f0, 2f0), φ=nothing, hZ=1f0 )
        cm_in  = CSG.flip(cm_out)

        pt_top    = CartesianPoint{Float32}(0, 1.5f0, 1f0)      # near top radius
        pt_bottom = CartesianPoint{Float32}(0, 0.5f0, -1f0)     # near bottom radius
        pt_mid    = CartesianPoint{Float32}(1f0, 0, 0f0)        # mid cone height

        normal_out = CSG.normal(cm_out, pt_mid)
        normal_in  = CSG.normal(cm_in, pt_mid)

        @test isa(normal_out, CartesianVector{Float32})
        @test isa(normal_in, CartesianVector{Float32})
        # inwards vs outwards should point opposite
        @test dot(normal_out, normal_in) ≈ -1f0 atol=1e-5
        
        n_arc = 8
        verts = CSG.vertices(cm_out, n_arc)
        @test isa(verts, Vector{CartesianPoint{Float32}})
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
        @test isa(pts, Vector{CartesianPoint{Float32}})
        @test length(pts) >= 2  # must have at least 2 points
        
        # Check that all points are within height bounds
        for pt in pts
            @test -cm_out.hZ - 1e-5 ≤ pt.z ≤ cm_out.hZ + 1e-5
        end
    end
    @testset "Ellipse" begin
        ell1 = @inferred CSG.Ellipse{Float32}(r = 1f0, φ=10f0)
        @inferred CSG.Ellipse{T}(r = (1,2))
        ell2 = @inferred CSG.Ellipse(r = 1f0, φ=10f0) 
        @test ell1 == ell2
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
        @test CSG.distance(pt_on, e) ≈ T(0)
        
        # Case 2: Point beyond b → distance = ||pt - b||
        pt_beyond = CartesianPoint{T}(2,1,0)
        @test CSG.distance(pt_beyond, e) ≈ norm(pt_beyond - b)
        
        # Case 3: Point before a → distance = ||pt - a||
        pt_before = CartesianPoint{T}(-1,1,0)
        @test CSG.distance(pt_before, e) ≈ norm(pt_before - a)
        
        # Case 4: Point off the segment but projection falls on segment
        pt_off = CartesianPoint{T}(0.5,2,0)
        v = b - a
        expected = norm(cross(pt_off - a, v)) / norm(v)
        @test CSG.distance(pt_off, e) ≈ expected
        
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
