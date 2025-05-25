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
        CSG.EllipsoidMantle{Float32}(r = 1f0)
        CSG.EllipsoidMantle() 
	  end
    @testset "EllipticalSurface" begin
        @inferred CSG.EllipticalSurface{Float32}(r = 1f0)
        @inferred CSG.EllipticalSurface()
        @inferred CSG.EllipticalSurface(φ=(1.,2.))    
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
    end
    @testset "Edge" begin
        edge1 = @inferred CSG.Edge{Float32}()
        edge2 = @inferred CSG.Edge(a = CartesianPoint{Float32}(0,0,0), b = CartesianPoint{Float16}(0,0,1)) 
        @test edge1 == edge2
    end
    @testset "Vector" begin
        cart = @inferred CSG.CartesianVector(x=2f0,z=1f0)
        @inferred CSG.CartesianVector{Float32}(x=2)
        cyl = @inferred CSG.CylindricalVector{Float32}(r=2.,z=1.)
        cyl2 = @inferred CSG.CylindricalVector(φ=3π)
        @test CartesianVector(cyl) == cart
        @test cart.x == Float32(2)
    end
    @testset "Point" begin
        cart = @inferred CSG.CartesianPoint(x=2f0,z=1f0)
        @inferred CSG.CartesianPoint{Float32}(x=2)
        cyl = @inferred CSG.CylindricalPoint{Float32}(r=2.,z=1.)
        cyl2 = @inferred CSG.CylindricalPoint(φ=3π)
        @test CartesianPoint(cyl) == cart
        @test cyl2.φ == T(π)
    end
end

@testset "Test geometry parsing" begin
    sim = Simulation{T}(SSD_examples[:InvertedCoax])
    c = sim.detector.contacts[2].geometry
    @test Geometry(T, Dictionary(c), default_units, no_translations) == c
end