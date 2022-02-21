# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).
using Test

# using CUDAKernels # Uncomment this line in order to run all tests on (CUDA) GPU
using SolidStateDetectors

using Unitful
using StaticArrays
using Tables, TypedTables

T = Float32
device_array_type = (@isdefined CUDAKernels) ? CUDAKernels.CUDA.CuArray : Array

@testset "Comparison to analytic solutions" begin
    include("comparison_to_analytic_solutions.jl")
end

@testset "SOR GPU Backend" begin
    include("SOR_GPU_Backend.jl")
end

@testset "ConstructiveSolidGeometry" begin
    include("ConstructiveSolidGeometry/CSG_IO.jl")
    include("ConstructiveSolidGeometry/CSG_primitives.jl")
end

T = Float32

@testset "Test real detectors" begin
    @testset "Simulate example detector: Inverted Coax" begin
        sim = Simulation{T}(SSD_examples[:InvertedCoax])
        simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 10e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].value[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
        nt = NamedTuple(sim)
        @test sim == Simulation(nt)
    end
    @testset "Simulate example detector: Inverted Coax (in cryostat)" begin
        sim = Simulation{T}(SSD_examples[:InvertedCoaxInCryostat])
        simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 10e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].value[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
    # @testset "Simulate example detector: Coax" begin
    #     sim = Simulation{T}(SSD_examples[:Coax])
    #     calculate_electric_potential!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
    #     calculate_electric_field!(sim)
    #     calculate_weighting_potential!(sim, 1, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
    #     calculate_weighting_potential!(sim, 2, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
    #     for i in 3:19
    #         calculate_weighting_potential!(sim, i, convergence_limit = 5e-6, refinement_limits = [0.2], verbose = false)
    #     end
    #     evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(30), 12e-3 )]))
    #     simulate!(evt, sim, Δt = 5e-10, max_nsteps = 10000)
    #     signalsum = T(0)
    #     for i in 1:length(evt.waveforms)
    #         signalsum += abs(ustrip(evt.waveforms[i].value[end]))
    #     end
    #     signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
    #     @info signalsum
    #     @test isapprox( signalsum, T(2), atol = 5e-3 )
    # end
    @testset "Simulate example detector: BEGe" begin
        sim = Simulation{T}(SSD_examples[:BEGe])
        calculate_electric_potential!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        calculate_electric_field!(sim)
        calculate_weighting_potential!(sim, 1, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        for i in 2:5
            calculate_weighting_potential!(sim, i, convergence_limit = 5e-6, device_array_type = device_array_type, refinement_limits = [0.2], verbose = false)
        end
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 20e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].value[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
        nt = NamedTuple(sim)
        @test sim == Simulation(nt)
    end
    @testset "Simulate example detector: HexagonalPrism" begin
        sim = Simulation{T}(SSD_examples[:Hexagon])
        simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event([CartesianPoint{T}(0, 5e-4, 1e-3)])
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].value[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
    @testset "Simulate example detector: CGD" begin
        sim = Simulation{T}(SSD_examples[:CGD])
        simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event([CartesianPoint{T}(0,2e-3,0)])
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].value[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
        nt = NamedTuple(sim)
        @test sim == Simulation(nt)
    end
    @testset "Simulate example detector: Spherical" begin
        sim = Simulation{T}(SSD_examples[:Spherical])
        simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event([CartesianPoint{T}(0,0,0)])
        simulate!(evt, sim)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].value[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end 
    @testset "Simulate example detector: Toroidal" begin
        sim = Simulation{T}(SSD_examples[:CoaxialTorus])
        SolidStateDetectors.apply_initial_state!(sim, ElectricPotential)
        simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1, 0.05, 0.02, 0.01], 
            max_tick_distance = 0.5u"mm", verbose = false)
        evt = Event([CartesianPoint{T}(0.0075,0,0)])
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].value[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
    @testset "Simulate example detector: SigGen PPC" begin
        sim = Simulation{T}(SSD_examples[:SigGen])
        simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 10e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].value[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
end

@testset "Diffusion and Self-Repulsion" begin
    sim = Simulation{T}(SSD_examples[:InvertedCoax])
    simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)

    pos = CartesianPoint{T}(0.02,0,0.05); Edep = 1u"eV"
    nbcc = NBodyChargeCloud(pos, Edep, 40, radius = T(0.0005), number_of_shells = 2)

    evt = Event(nbcc)
    simulate!(evt, sim, self_repulsion = true, diffusion = true)
    signalsum = T(0)
    for i in 1:length(evt.waveforms)
        signalsum += abs(ustrip(evt.waveforms[i].value[end]))
    end
    signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
    @info signalsum
    @test isapprox( signalsum, T(2), atol = 5e-3 )
end

@testset "Table Simulation" begin 
    sim = Simulation{T}(SSD_examples[:InvertedCoax])
    simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)

    evt_table = Table(
        evtno = Int32[1], 
        detno = Int32[1],
        thit = [T[0] * u"s"],
        edep = [T[1] * u"eV"],
        pos = [[SVector{3, T}.(0.01, 0.01, 0.01) * u"m"]]
    )
    contact_charge_signals = simulate_waveforms(      
        evt_table,
        sim,
        max_nsteps = 4000, 
        Δt = 1u"ns", 
        number_of_carriers = 20,
        number_of_shells = 2,
        verbose = false);
    signalsum = T(0)
    for i in 1:length(contact_charge_signals.waveform)
        signalsum += abs(ustrip(contact_charge_signals.waveform[i].value[end]))
    end
    signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
    @test isapprox( signalsum, T(2), atol = 5e-3 )
end

@testset "ADLChargeDriftModel" begin
    include("ADLChargeDriftModel.jl")
end

@testset "IO" begin
    sim = Simulation(SSD_examples[:InvertedCoax])
    
    calculate_electric_potential!(sim, verbose = false, device_array_type = device_array_type)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)
    
    calculate_electric_field!(sim)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    calculate_weighting_potential!(sim, 1, verbose = false, device_array_type = device_array_type)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    for i in findall(ismissing.(sim.weighting_potentials)) calculate_weighting_potential!(sim, i, verbose = false, device_array_type = device_array_type) end
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)
end 


