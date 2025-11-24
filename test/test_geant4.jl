# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors
using ArraysOfArrays
using RadiationDetectorSignals
using StaticArrays
using TypedTables
using Unitful

import Geant4

T = Float32

@testset "Parse materials in Geant4 extension" begin
    parse_material = Base.get_extension(SolidStateDetectors, :SolidStateDetectorsGeant4Ext).parse_material
    for (m,material) in SolidStateDetectors.material_properties
        @testset "$(m)" begin
            @test_nowarn parse_material(material.name)
        end
    end
end

@testset "Define particle sources" begin

    # Isotropic emission when no direction is passed
    m1 = MonoenergeticSource("gamma", 2.615u"MeV", CartesianPoint(0.04, 0, 0.05))
    @test m1.direction == CartesianVector(0,0,1)
    @test m1.opening_angle == 180u"°"

    # Directed ray emission when no opening_angle is passed
    m2 = MonoenergeticSource("gamma", 2.615u"MeV", CartesianPoint(0.04, 0, 0.05), CartesianVector(-1, 0, 0))
    @test m2.direction == CartesianVector(-1, 0, 0)
    @test m2.opening_angle == 0u"°"

    # Directed cone emission when opening_angle is passed
    m3 = MonoenergeticSource("gamma", 2.615u"MeV", CartesianPoint(0.04, 0, 0.05), CartesianVector(-1, 0, 0), 10u"°")
    @test m3.direction == CartesianVector(-1, 0, 0)
    @test m3.opening_angle == 10u"°"

    # Isotropic emission when no direction is passed
    i1 = IsotopeSource(82, 212, 0, 0, CartesianPoint(0.04, 0, 0.05))
    @test i1.direction == CartesianVector(0,0,1)
    @test i1.opening_angle == 180u"°"

    # Directed ray emission when no opening_angle is passed
    i2 = IsotopeSource(82, 212, 0, 0, CartesianPoint(0.04, 0, 0.05), CartesianVector(-1, 0, 0))
    @test i2.direction == CartesianVector(-1, 0, 0)
    @test i2.opening_angle == 0u"°"

    # Directed cone emission when opening_angle is passed
    i3 = IsotopeSource(82, 212, 0, 0, CartesianPoint(0.04, 0, 0.05), CartesianVector(-1, 0, 0), 10u"°")
    @test i3.direction == CartesianVector(-1, 0, 0)
    @test i3.opening_angle == 10u"°"

end

@testset "Run Geant4 simulations" begin

    # Create G4JLApplication
    sim = @test_nowarn Simulation{T}(SSD_examples[:InvertedCoaxInCryostat])
    source = MonoenergeticSource("gamma", 2.615u"MeV", CartesianPoint(0.04, 0, 0.05), CartesianVector(-1,0,0))
    app = Geant4.G4JLApplication(sim, source)

    # Simulate 100 events
    evts = run_geant4_simulation(app, 100)
    @test evts isa Table
    @test SolidStateDetectors.is_detector_hits_table(evts)
    @test length(evts) == 100

    # Add fano noise
    material = SolidStateDetectors.material_properties[:HPGe]
    evts_fano = add_fano_noise(evts, material.E_ionisation, material.f_fano)
    @test evts_fano isa Table
    @test SolidStateDetectors.is_detector_hits_table(evts_fano)
    @test length(evts_fano) == 100

    # Cluster events by radius
    clustered_evts = SolidStateDetectors.cluster_detector_hits(evts, 10u"µm")
    @test length(clustered_evts) == length(evts)
    @test length(flatview(clustered_evts.pos)) <= length(flatview(evts.pos))
    @test eltype(first(clustered_evts.pos)) <: StaticVector{3}

    # Generate waveforms
    simulate!(sim, refinement_limits = [0.2,0.1,0.05,0.03,0.02])
    wf = simulate_waveforms(evts, sim, Δt = 1u"ns", max_nsteps = 2000)
    @test wf isa Table
    @test :waveform in columnnames(wf)
    @test length(wf) == length(evts) * sum(.!ismissing.(sim.weighting_potentials))

    # Use the table with added Fano noise
    wf = simulate_waveforms(evts_fano, sim, Δt = 1u"ns", max_nsteps = 2000)
    @test wf isa Table
    @test :waveform in columnnames(wf)
    @test length(wf) == length(evts_fano) * sum(.!ismissing.(sim.weighting_potentials))
end