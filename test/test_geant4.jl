# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors
using ArraysOfArrays
using RadiationDetectorSignals
using Rotations
using StaticArrays
using TypedTables
using Unitful
using LegendHDF5IO

import Geant4

T = Float32

@testset "Parse materials in Geant4 extension" begin
    parse_material = Base.get_extension(SolidStateDetectors, :SolidStateDetectorsGeant4Ext).parse_material
    for (m,material) in SolidStateDetectors.material_properties
        @testset "$(m)" begin
            @test @test_nowarn parse_material(material.name) isa AbstractString
        end
    end
    @test_throws ArgumentError parse_material("Unknown material")
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
    sim = Simulation{T}(SSD_examples[:InvertedCoaxInCryostat])
    source = MonoenergeticSource("gamma", 2.615u"MeV", CartesianPoint(0.04, 0, 0.05), CartesianVector(-1,0,0))
    app = Geant4.G4JLApplication(sim, source, verbose = false)

    # Simulate Geant4 events
    N = 100
    evts = run_geant4_simulation(app, N)
    @test evts isa Table
    @test SolidStateDetectors.is_detector_hits_table(evts)
    @test SolidStateDetectors.get_detector_hits_table(evts) == evts
    @test length(evts) == N
    @test all(length.(values(columns(evts))) .== N)

    # Add fano noise
    material = SolidStateDetectors.material_properties[:HPGe]
    evts_fano = add_fano_noise(evts, material.E_ionisation, material.f_fano)
    @test evts_fano isa Table
    @test SolidStateDetectors.is_detector_hits_table(evts_fano)
    @test length(evts_fano) == N

    # Cluster events by radius
    clustered_evts = SolidStateDetectors.cluster_detector_hits(evts, 10u"µm")
    @test length(clustered_evts) >= length(evts)
    @test length(flatview(clustered_evts.pos)) <= length(flatview(evts.pos))
    column_lengths = unique(length.(values(columns(clustered_evts))))
    @test length(column_lengths) == 1
    @test only(column_lengths) >= N
    @test eltype(first(clustered_evts.pos)) <: CartesianPoint

    # Cluster events by radius, but split in time
    clustered_evts = SolidStateDetectors.cluster_detector_hits(evts, 10u"µm", 1u"ns")
    @test length(evts) <= length(clustered_evts) <= length(flatview(evts.pos))
    @test length(flatview(clustered_evts.pos)) <= length(flatview(evts.pos))
    column_lengths = unique(length.(values(columns(clustered_evts))))
    @test length(column_lengths) == 1
    @test only(column_lengths) >= N
    @test eltype(first(clustered_evts.pos)) <: CartesianPoint

    # Try the same method using StaticVectors as eltype of pos
    evts_static = Table(evts; pos = VectorOfVectors(broadcast.(p -> SVector{3}(p.x, p.y, p.z), evts.pos)))
    clustered_evts_static = SolidStateDetectors.cluster_detector_hits(evts_static, 10u"µm", 1u"ns")
    @test length(clustered_evts_static) == length(evts_static)
    @test length(flatview(clustered_evts_static.pos)) <= length(flatview(evts_static.pos))
    column_lengths = unique(length.(values(columns(clustered_evts_static))))
    @test length(column_lengths) == 1
    @test only(column_lengths) >= N
    @test eltype(first(clustered_evts_static.pos)) <: StaticVector{3}

    # Temporary HDF5 file
    mktemp() do tmpfile, io
        LegendHDF5IO.lh5open(tmpfile, "w") do h5
            LegendHDF5IO.writedata(h5, "geant4", evts)
        end
        @test isfile(tmpfile)
        @test evts == LegendHDF5IO.lh5open(tmpfile) do h5
            LegendHDF5IO.readdata(h5, "geant4")
        end
    end

    # Cluster aggressively and check idempotence
    combine_evts = SolidStateDetectors.cluster_detector_hits(evts, 1u"m")
    @test length(combine_evts) == length(evts)
    @test SolidStateDetectors.cluster_detector_hits(combine_evts, 1u"m") == combine_evts
    
    # Generate waveforms
    timed_simulate!(sim, refinement_limits = [0.2,0.1,0.05,0.03,0.02])
    evt = Event(evts[1], T)
    @test evt isa Event{T}
    timed_simulate!(evt, sim, Δt = 1u"ns", max_nsteps = 2000)

    # Test generating electron and hole contribution separately
    (; electron_contribution, hole_contribution) = get_electron_and_hole_contribution(evt, sim, 1)
    @test electron_contribution isa RDWaveform
    @test hole_contribution isa RDWaveform

    # Use the raw table
    wf = timed_simulate_waveforms(evts, sim, Δt = 1u"ns", max_nsteps = 2000)
    @test wf isa Table
    @test SolidStateDetectors.get_detector_hits_table(wf[findall(wf.chnid .== 1)]) == evts
    @test :waveform in columnnames(wf)
    @test length(wf) == length(evts) * sum(.!ismissing.(sim.weighting_potentials))
    @test wf.waveform[1] ≈ evt.waveforms[1]

    # Generate waveforms from the clustered events
    wf_clustered = timed_simulate_waveforms(clustered_evts, sim, Δt = 1u"ns", max_nsteps = 2000)
    @test wf_clustered isa Table
    @test :waveform in columnnames(wf_clustered)
    @test length(wf_clustered) == length(clustered_evts) * sum(.!ismissing.(sim.weighting_potentials))

    # Use the table with added Fano noise
    wf_fano = timed_simulate_waveforms(evts_fano, sim, Δt = 1u"ns", max_nsteps = 2000)
    @test wf_fano isa Table
    @test :waveform in columnnames(wf_fano)
    @test length(wf_fano) == length(evts_fano) * sum(.!ismissing.(sim.weighting_potentials))

    # Try the same method using StaticVectors as eltype of pos (with units of mm)
    evts_static = Table(evts; pos = VectorOfVectors(broadcast.(p -> SVector{3}(p.x, p.y, p.z) * 1000u"mm", SolidStateDetectors.to_internal_units(evts.pos))))
    clustered_evts_static = SolidStateDetectors.cluster_detector_hits(evts_static, 10u"µm", 1u"ns")
    @test length(clustered_evts_static) == length(evts_static)
    @test length(flatview(clustered_evts_static.pos)) <= length(flatview(evts_static.pos))
    @test eltype(first(clustered_evts_static.pos)) <: StaticVector{3}

    # Check clustering with cluster radius being unitless or with the wrong unit
    clustered_evts_static_unitless = SolidStateDetectors.cluster_detector_hits(evts_static, 1e-5, 1u"ns")
    @test clustered_evts_static == clustered_evts_static_unitless
    clustered_evts_static_unitless = SolidStateDetectors.cluster_detector_hits(evts_static, 10u"µm", 1e-9)
    @test clustered_evts_static == clustered_evts_static_unitless
    @test_throws Exception SolidStateDetectors.cluster_detector_hits(evts_static, 2u"s", 1u"ns")
    @test_throws Exception SolidStateDetectors.cluster_detector_hits(evts_static, 10u"µm", 2u"mm")

    # Generate waveforms using StaticVectors as eltype of pos
    wf_static = timed_simulate_waveforms(evts_static, sim, Δt = 1u"ns", max_nsteps = 2000)
    @test wf_static isa Table
    @test :waveform in columnnames(wf_static)
    @test length(wf_static) == length(evts_static) * sum(.!ismissing.(sim.weighting_potentials))

    # Result should not depend whether input positions are SVector or CartesianPoint
    @test wf_static.waveform ≈ wf.waveform

    # test table_utils.jl
    @testset "Test table utils" begin
        evts_split = SolidStateDetectors.split_table_by_each_charge_deposition(evts)
        @test length(evts_split) == length(flatview(evts.pos))
            
        v = CartesianVector(one(T)*u"m", zero(T)*u"m", zero(T)*u"m")
        v_static = SVector(T(1000)*u"mm", zero(T)*u"mm", zero(T)*u"mm")
        translated = SolidStateDetectors.translate_event_positions(evts, v)
        @test all(p -> p isa eltype(first(evts.pos)), flatview(translated.pos))
        @test map(p -> p + v, flatview(evts.pos)) == flatview(translated.pos)
        translated_static = SolidStateDetectors.translate_event_positions(evts_static, v_static)
        @test all(p -> p isa eltype(first(evts_static.pos)), flatview(translated_static.pos))
        @test map(p -> p + v_static, flatview(evts_static.pos)) == flatview(translated_static.pos)

        # inter-compatibility
        t1 = SolidStateDetectors.translate_event_positions(evts_static, v)
        @test all(p -> p isa eltype(first(evts_static.pos)), flatview(t1.pos))
        @test all(isapprox.(flatview(translated_static.pos), flatview(t1.pos)))
        t2 = SolidStateDetectors.translate_event_positions(evts, v_static)
        @test all(p -> p isa eltype(first(evts.pos)), flatview(t2.pos))
        @test all(isapprox.(flatview(translated.pos), flatview(t2.pos)))

        r = RotZ(π)
        rotated = SolidStateDetectors.rotate_event_positions(evts, r)
        @test all(p -> p isa eltype(first(evts.pos)), flatview(rotated.pos))
        @test map(p -> eltype(first(evts.pos))(cartesian_zero + r * (p - cartesian_zero)) , flatview(evts.pos)) == flatview(rotated.pos)
        rotated_static = SolidStateDetectors.rotate_event_positions(evts_static, r)
        @test all(p -> p isa eltype(first(evts_static.pos)), flatview(rotated_static.pos))
        @test map(p -> eltype(first(evts_static.pos))(r * p), flatview(evts_static.pos)) == flatview(rotated_static.pos)
    end

    # Test IsotopeSource (Cs-137)
    isource = IsotopeSource(55, 137, 0, 0, CartesianPoint(0.05,0,0.05))
    iapp = Geant4.G4JLApplication(sim, isource, verbose = false)
    iapp.sdetectors = app.sdetectors # workaround: copy sdetectors from previous G4JLApplication to avoid undefined reference
    ievts = run_geant4_simulation(iapp, N)

    @test ievts isa Table
    @test SolidStateDetectors.is_detector_hits_table(ievts)
    @test length(ievts) == N

end

@testset "Construct G4JLDetectors from config files" begin
    for (key, sim_file) in SSD_examples
        if endswith(sim_file, ".yaml")
            @testset "Read from SSD config file: $(key)" begin
                fn = tempname() * ".gdml"
                @test Geant4.G4JLDetector(sim_file, fn, save_gdml = true, verbose = false) isa Geant4.G4JLDetector
                @test isfile(fn)
                @test @test_nowarn Geant4.G4JLDetector(fn) isa Geant4.G4JLDetector
                rm(fn)
            end
        end
    end
end