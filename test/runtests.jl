# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).
using Test

using SolidStateDetectors

using Unitful

T = Float32

@testset "Comparison to analytic solutions" begin
    include("comparison_to_analytic_solutions.jl")
end

@testset "Test real detectors" begin
    @testset "Simulate example detector: Inverted Coax" begin
        sim = Simulation{T}(SSD_examples[:InvertedCoax])
        simulate!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 10e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
        nt = NamedTuple(sim)
        @test sim == Simulation(nt)
    end
    @testset "Simulate example detector: Inverted Coax (in cryostat)" begin
        sim = Simulation{T}(SSD_examples[:InvertedCoaxInCryostat])
        simulate!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 10e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
    @testset "Simulate example detector: Coax" begin
        sim = Simulation{T}(SSD_examples[:Coax])
        calculate_electric_potential!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
        calculate_electric_field!(sim)
        calculate_drift_fields!(sim)
        calculate_weighting_potential!(sim, 1, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
        calculate_weighting_potential!(sim, 2, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
        for i in 3:19
            calculate_weighting_potential!(sim, i, convergence_limit = 5e-6, refinement_limits = [0.2], verbose = false)
        end
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(30), 12e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
    @testset "Simulate example detector: BEGe" begin
        sim = Simulation{T}(SSD_examples[:BEGe])
        calculate_electric_potential!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
        calculate_electric_field!(sim)
        calculate_drift_fields!(sim)
        calculate_weighting_potential!(sim, 1, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
        for i in 2:5
            calculate_weighting_potential!(sim, i, convergence_limit = 5e-6, refinement_limits = [0.2], verbose = false)
        end
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 20e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
        nt = NamedTuple(sim)
        @test sim == Simulation(nt)
    end
    @testset "Simulate example detector: HexagonalPrism" begin
        sim = Simulation{T}(SSD_examples[:Hexagon])
        simulate!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event([CartesianPoint{T}(0, 5e-4, 1e-3)])
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
    @testset "Simulate example detector: CGD" begin
        sim = Simulation{T}(SSD_examples[:CGD])
        simulate!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event([CartesianPoint{T}(0,2e-3,0)])
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
        nt = NamedTuple(sim)
        @test sim == Simulation(nt)
    end
    @testset "Simulate example detector: Spherical" begin
        sim = Simulation{T}(SSD_examples[:Spherical])
        simulate!(sim, convergence_limit = 1e-5, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event([CartesianPoint{T}(0,0,0)])
        simulate!(evt, sim)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end 
    @testset "Simulate example detector: Toroidal" begin
        sim = Simulation{T}(SSD_examples[:CoaxialTorus])
        SolidStateDetectors.apply_initial_state!(sim, ElectricPotential)
        simulate!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1, 0.05], 
            max_tick_distance = 0.5u"mm", verbose = false)
        evt = Event([CartesianPoint{T}(0.01,0,0.003)])
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
    @testset "Simulate example detector: SigGen PPC" begin
        sim = Simulation{T}(SSD_examples[:SigGen])
        simulate!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 10e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
end

@testset "ADLChargeDriftModel" begin
    include("ADLChargeDriftModel.jl")
end

@testset "IO" begin
    sim = Simulation(SSD_examples[:InvertedCoax])
    
    calculate_electric_potential!(sim, verbose = false)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)
    
    calculate_electric_field!(sim)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)
    
    calculate_drift_fields!(sim)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    calculate_weighting_potential!(sim, 1, verbose = false)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    for i in findall(ismissing.(sim.weighting_potentials)) calculate_weighting_potential!(sim, i, verbose = false) end
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)
end 

include("ConstructiveSolidGeometry/CSG_IO.jl")
