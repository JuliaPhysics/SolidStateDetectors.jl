using Test
using SolidStateDetectors

import SolidStateDetectors.ConstructiveSolidGeometry as CSG
import SolidStateDetectors.ConstructiveSolidGeometry: Geometry

import LegendHDF5IO
import Unitful: @u_str
import ArraysOfArrays: VectorOfVectors

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

@testset "Test LegendHDF5IO with AbstractCoordinatePoint" begin
    x = rand(T, 3)
    x_unit = CSG.from_internal_units(x, u"mm")
    pt = CartesianPoint(x...)
    pt_unit = CartesianPoint(x_unit...)
    pt_cyl = CylindricalPoint(x...)
    pt_cyl_unit = CylindricalPoint(CSG.from_internal_units.(x, [u"mm", u"°", u"mm"])...)

    @testset begin
        mktemp() do tmpfile, io
            for x in (x, x_unit, pt, pt_unit, pt_cyl)    
                LegendHDF5IO.lh5open(tmpfile, "w") do h
                    
                    @testset "I/O for AbstractCoordinatePoints" begin
                        # write to file using writedata and LH5Array
                        @test_nowarn LegendHDF5IO.writedata(h, "x1", x)
                        @test_nowarn h["x2"] = x

                        # read using readddata
                        x1 = @test_nowarn LegendHDF5IO.readdata(h, "x1")
                        @test x1 isa typeof(x)
                        @test x1 == x

                        # read using LH5Array
                        x2 = @test_nowarn h["x2"]
                        @test x2 == x
                        
                        # cross-compatibility between the two methods
                        x1_ = @test_nowarn h["x1"]
                        x2_ = @test_nowarn LegendHDF5IO.readdata(h, "x2")
                        @test x1_ == x2_ == x
                    end
                    
                    @testset "I/O for Vector of AbstractCoordinatePoints" begin
                        
                        y = x isa Vector ? VectorOfVectors([x]) : [x]

                        # write to file using writedata and LH5Array
                        @test_nowarn LegendHDF5IO.writedata(h, "y1", y)
                        @test_nowarn h["y2"] = y

                        # read using readddata
                        y1 = @test_nowarn LegendHDF5IO.readdata(h, "y1")
                        @test y1 isa typeof(y)
                        @test y == y1

                        # read using LH5Array
                        @test y == h["y2"]
                        
                        # cross-compatibility between the two methods
                        y1_ = @test_nowarn h["y1"]
                        y2_ = @test_nowarn LegendHDF5IO.readdata(h, "y2")
                        @test y1_ == y2_ == y
                    end
                end
            end
        end
    end
end