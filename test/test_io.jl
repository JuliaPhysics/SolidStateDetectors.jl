# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors

using ArraysOfArrays
using LegendHDF5IO
using OrderedCollections
using StaticArrays
using TypedTables
using Unitful

# include("test_utils.jl")

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

    @test_logs (:info,) timed_calculate_electric_field!(sim, n_points_in_φ = 15)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    timed_calculate_weighting_potential!(sim, 1, verbose = false)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    for i in findall(ismissing.(sim.weighting_potentials)) timed_calculate_weighting_potential!(sim, i, verbose = false) end
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    for (Potential, spot) in (
            (ElectricPotential, sim.electric_potential), 
            (ElectricField, sim.electric_field), 
            (ImpurityScale, sim.imp_scale), 
            (PointTypes, sim.point_types), 
            (EffectiveChargeDensity, sim.q_eff_imp), 
            (EffectiveChargeDensity, sim.q_eff_fix), 
            (WeightingPotential, sim.weighting_potentials[1]), 
            (DielectricDistribution, sim.ϵ_r)
        )

        @test length(spot) isa Integer
        @test spot[1] isa eltype(spot.data)
        @test spot[:z] isa SolidStateDetectors.DiscreteAxis
        # test the implementation of Base.convert
        spot_nt::NamedTuple = spot
        spot2::Potential = spot_nt
        @test spot2 == spot
    end
end

@testset "Test LegendHDF5IO" begin
    fn = tempname() * ".lh5"
    @test_nowarn ssd_write(fn, sim)
    @test isfile(fn)
    @test_logs (:warn, "Destination `$(fn)` already exists. Overwriting...") ssd_write(fn, sim)
    @test sim == ssd_read(fn, Simulation)
    rm(fn)
end

@timed_testset "Table Simulation" begin 

    T = Float32 
    sim = Simulation{T}(SSD_examples[:InvertedCoax])
    timed_simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)

    # test I/O for both SVector (with and without units) and CartesianPoint
    for pos in ([SVector{3, T}.(10, 10, 10) * u"mm"], [SVector{3, T}.(0.01, 0.01, 0.01)], [CartesianPoint{T}(0.01, 0.01, 0.01)])
        
        evt_table = Table(
            evtno = Int32[1], 
            detno = Int32[1],
            thit = VectorOfVectors([T[0] * u"s"]),
            edep = VectorOfVectors([T[1] * u"eV"]),
            pos = VectorOfVectors([pos])
        )

        contact_charge_signals = timed_simulate_waveforms(      
                evt_table,
                sim,
                max_nsteps = 4000, 
                Δt = 1u"ns", 
                number_of_carriers = 20,
                number_of_shells = 2,
                geometry_check = true,
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
                geometry_check = true,
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
end
