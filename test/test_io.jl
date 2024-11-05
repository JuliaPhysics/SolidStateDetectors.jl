# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors

using ArraysOfArrays
using LegendHDF5IO
using OrderedCollections
using StaticArrays
using TypedTables
using Unitful

include("test_utils.jl")

@testset "Read in different config file formats" begin

    # Create JSON file from YAML file
    filename = SSD_examples[:InvertedCoax]
    @test_nowarn SolidStateDetectors.yaml2json(filename)
    
    # Test read in using YAML, JSON and SigGen using Dict, OrderedDict and OrderedDict{String, Any}
    for f in (filename, replace(filename, ".yaml" => ".json"), SSD_examples[:SigGen])
        @test isfile(f)
        @testset "$(split(f, ".")[end])" begin
            @test_nowarn SolidStateDetectors.parse_config_file(f)
            @test_nowarn SolidStateDetectors.parse_config_file(f, dicttype = OrderedDict)
            @test_nowarn SolidStateDetectors.parse_config_file(f, dicttype = OrderedDict{String, Any})
        end
    end
end

sim = Simulation(SSD_examples[:InvertedCoax])

@testset "Conversion Simulation <--> NamedTuple" begin
    timed_calculate_electric_potential!(sim, verbose = false)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    timed_calculate_electric_field!(sim)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    timed_calculate_weighting_potential!(sim, 1, verbose = false)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    for i in findall(ismissing.(sim.weighting_potentials)) timed_calculate_weighting_potential!(sim, i, verbose = false) end
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)
end

@testset "Test LegendHDF5IO" begin
    @test_nowarn ssd_write("legendhdf5io_test.lh5", sim)
    @test isfile("legendhdf5io_test.lh5")
    @test_logs (:warn, "Destination `legendhdf5io_test.lh5` already exists. Overwriting...") ssd_write("legendhdf5io_test.lh5", sim)
    @test sim == ssd_read("legendhdf5io_test.lh5", Simulation)
    rm("legendhdf5io_test.lh5")
end

@timed_testset "Table Simulation" begin 

    T = Float32 
    sim = Simulation{T}(SSD_examples[:InvertedCoax])
    timed_simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)

    evt_table = Table(
        evtno = Int32[1], 
        detno = Int32[1],
        thit = VectorOfVectors([T[0] * u"s"]),
        edep = VectorOfVectors([T[1] * u"eV"]),
        pos = VectorOfVectors([[SVector{3, T}.(0.01, 0.01, 0.01) * u"m"]])
    )

    contact_charge_signals = timed_simulate_waveforms(      
            evt_table,
            sim,
            max_nsteps = 4000, 
            Δt = 1u"ns", 
            number_of_carriers = 20,
            number_of_shells = 2,
            verbose = false);
    
    @timed_testset "Regular simulate_waveforms" begin
        signalsum = T(0)
        for i in 1:length(contact_charge_signals.waveform)
            signalsum += abs(ustrip(contact_charge_signals.waveform[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end

    @timed_testset "LegendHDF5IO simulate_waveforms" begin
        timed_simulate_waveforms(      
            evt_table,
            sim,
            ".",
            chunk_n_physics_events = 1,
            max_nsteps = 4000, 
            Δt = 1u"ns", 
            number_of_carriers = 20,
            number_of_shells = 2,
            verbose = false
        )
        @info isfile("generated_waveforms_evts_1-1.h5")
        @test contact_charge_signals == LegendHDF5IO.lh5open("generated_waveforms_evts_1-1.h5") do h5f
            @test haskey(h5f, "generated_waveforms")
            LegendHDF5IO.readdata(h5f.data_store, "generated_waveforms")
        end
        rm("generated_waveforms_evts_1-1.h5")
    end
end