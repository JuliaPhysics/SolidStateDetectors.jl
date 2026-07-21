# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using SolidStateDetectors
using Test
using Unitful

import SolidStateDetectors: ConfigFileError, parse_config_file, construct_units

T = Float32

@testset "Invalid grid coordinates raise ConfigFileError" begin
    cfg = parse_config_file(SSD_examples[:InvertedCoax])
    cfg["grid"]["coordinates"] = "spherical"
    @test_throws ConfigFileError Simulation{T}(cfg)
end

@testset "Unknown charge trapping model raises ConfigFileError" begin
    cfg = parse_config_file(SSD_examples[:InvertedCoax])
    cfg["detectors"][1]["semiconductor"]["charge_trapping_model"] = Dict("model" => "Bogss")
    @test_throws ConfigFileError Simulation{T}(cfg)
end

@testset "Doped-contact auto-detection without doped_contact_id" begin
    cfg = parse_config_file(SSD_examples[:IVCIlayer])
    delete!(cfg["detectors"][1]["semiconductor"]["impurity_density"], "doped_contact_id")
    det = SolidStateDetector{T}(cfg, construct_units(cfg))
    @test det isa SolidStateDetector
end

@testset "Celsius temperature units parse consistently" begin
    cfg = parse_config_file(SSD_examples[:IVCIlayer])
    cfg["units"]["temperature"] = "Celsius"
    temp_C = -183.15  # 90 K
    cfg["detectors"][1]["semiconductor"]["temperature"] = temp_C
    cfg["detectors"][1]["semiconductor"]["charge_drift_model"]["temperature"] = temp_C
    sim = Simulation{T}(cfg)
    @test sim.detector.semiconductor.temperature ≈ T(90) atol = T(0.01)
    @test sim.detector.semiconductor.charge_drift_model.temperature ≈ T(90) atol = T(0.01)
end

@testset "2D and phi-reduced grids require the matching detector symmetry" begin
    # cuboid detector on a cylindrical grid: 2D must be rejected
    cfg = parse_config_file(SSD_examples[:CGD_CylGrid])
    cfg["grid"]["axes"]["phi"] = Dict("from" => 0, "to" => 0, "boundaries" => "periodic")
    sim2d = Simulation{T}(cfg)
    @test_throws ConfigFileError Grid(sim2d)

    # ... and a periodic 120° wedge of it as well (no 3-fold symmetry)
    cfg["grid"]["axes"]["phi"] = Dict("from" => 0, "to" => 120, "boundaries" => "periodic")
    simw = Simulation{T}(cfg)
    @test_throws ConfigFileError Grid(simw)

    # φ-symmetric detector in 2D stays fine
    @test Grid(Simulation{T}(SSD_examples[:InvertedCoax])) isa SolidStateDetectors.CylindricalGrid{T}

    # the segmented BEGe example (3-fold symmetric, reduced φ range) stays fine
    @test Grid(Simulation{T}(SSD_examples[:BEGe])) isa SolidStateDetectors.CylindricalGrid{T}
end

@testset "2D symmetry guard rejects discretely-but-not-continuously symmetric detectors" begin
    # 8-pointed star (two boxes rotated by 45°): C8-symmetric but not axisymmetric,
    # so it must not pass as a 2D detector despite being invariant under all π/4 shifts
    cfg = parse_config_file(SSD_examples[:CGD_CylGrid])
    cfg["grid"]["axes"]["phi"] = Dict("from" => 0, "to" => 0, "boundaries" => "periodic")
    cfg["detectors"][1]["semiconductor"]["geometry"] = Dict(
        "union" => [
            Dict("box" => Dict("widths" => [10, 10, 10])),
            Dict("box" => Dict("widths" => [10, 10, 10], "rotate" => Dict("Z" => 45))),
        ])
    cfg["detectors"][1]["contacts"][1]["geometry"] = Dict(
        "tube" => Dict("r" => 3, "h" => 0.5, "origin" => Dict("z" => 5.25)))
    cfg["detectors"][1]["contacts"][2]["geometry"] = Dict(
        "tube" => Dict("r" => 3, "h" => 0.5, "origin" => Dict("z" => -5.25)))
    delete!(cfg, "surroundings")
    @test_throws ConfigFileError Grid(Simulation{T}(cfg))
end

# axisymmetric geometry with a φ-dependent impurity density: also not valid in 2D
struct PhiDependentImpurityDensity{T} <: SolidStateDetectors.AbstractImpurityDensity{T} end
function SolidStateDetectors.get_impurity_density(
        ::PhiDependentImpurityDensity{T}, pt::SolidStateDetectors.AbstractCoordinatePoint{T}) where {T}
    cpt = SolidStateDetectors.CylindricalPoint(pt)
    T(1e16) * (1 + cpt.r * cos(cpt.φ))
end

@testset "2D symmetry guard rejects φ-dependent impurity densities" begin
    sim = Simulation{T}(SSD_examples[:InvertedCoax])
    sim.detector = SolidStateDetector(sim.detector, PhiDependentImpurityDensity{T}())
    @test_throws ConfigFileError Grid(sim)
    @test Grid(sim; check_φ_symmetry = false) isa SolidStateDetectors.CylindricalGrid{T}
end
