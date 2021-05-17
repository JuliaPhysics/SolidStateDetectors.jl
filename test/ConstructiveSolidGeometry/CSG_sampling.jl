using Test

using SolidStateDetectors.ConstructiveSolidGeometry: sample,
    ConalPlane, ConeMantle, CylindricalAnnulus, 
    RectangleX, RectangleY, RectangleZ, 
    RegularHexagon, HexagonalPrismMantle, RegularPolygon, RegularPrismMantle,
    SphereMantle, ToroidalAnnulus, TorusMantle,
    Tube, Cone, Torus, Box, Sphere, HexagonalPrism

@testset "Test sampling" begin

    @testset "VolumePrimitives" begin
        
        volume_primitives = [Tube, Cone, Torus, Box, Sphere, HexagonalPrism]
        
        for P in volume_primitives
            @testset "$P" begin
                p = P(); @test all(broadcast(pt -> pt in p, sample(p)))
                p = P(); @test all(broadcast(pt -> pt in p, sample(p, (10,10,10))))
            end
        end
    end
    
    @testset "SurfacePrimitives" begin
        
        surface_primitives = [ConalPlane, CylindricalAnnulus,
                    #ConeMantle, RegularHexagon, HexagonalPrismMantle,
                    RectangleX, RectangleY, RectangleZ, 
                    SphereMantle, ToroidalAnnulus, TorusMantle]
        
        @testset "Cartesian Grid" begin
            cart = (x = collect(-2:0.01:2), y = collect(-2:0.01:2),   z = collect(-2:0.01:2))
            for P in surface_primitives
                @testset "$P" begin
                    p = P(); @test all(broadcast(pt -> pt in p, sample(p, cart)))
                end
            end
        end
        
        surface_primitives = [ConalPlane, CylindricalAnnulus,
                    ConeMantle, RegularHexagon, HexagonalPrismMantle,
                    #RectangleX, RectangleY, RectangleZ, 
                    SphereMantle, ToroidalAnnulus, TorusMantle]
        
        @testset "Cylindrical Grid" begin
            cylt = (r = collect(0:0.01:2),  φ = collect(0:pi/8:2pi), z = collect(-2:0.01:2))
            for P in surface_primitives
                @testset "$P" begin
                    p = P(); @test all(broadcast(pt -> pt in p, sample(p, cylt)))
                end
            end
        end
    end
end