# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).
using Test

using SolidStateDetectors
using SolidStateDetectors: SSDFloat

using Unitful

T = Float32

@testset "Comparison to analytic solutions" begin
    include("comparison_to_analytic_solutions.jl")
end

@testset "Test real detectors" begin
    @testset "Simulate example detector: Inverted Coax" begin
        sim = Simulation{T}(SSD_examples[:InvertedCoax])
        simulate!(sim, max_refinements = 1, verbose = true)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 10e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-4 )
    end
    @testset "Simulate example detector: Inverted Coax (in cryostat)" begin
        sim = Simulation{T}(SSD_examples[:InvertedCoaxInCryostat])
        simulate!(sim, max_refinements = 1, verbose = true)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 10e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-4 )
    end
    @testset "Simulate example detector: Coax" begin
        sim = Simulation{T}(SSD_examples[:Coax])
        simulate!(sim, max_refinements = 0, verbose = true)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(30), 12e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 0.02 )
    end
    @testset "Simulate example detector: BEGe" begin
        sim = Simulation{T}(SSD_examples[:BEGe])
        simulate!(sim, max_refinements = 1, verbose = true)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 20e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 3e-2 )
    end
    @testset "Simulate example detector: HexagonalPrism" begin
        sim = Simulation{T}(SSD_examples[:Hexagon])
        simulate!(sim, max_refinements = 0, verbose = true)
        evt = Event([CartesianPoint{T}(0, 5e-4, 1e-3)])
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-4 )
    end
    @testset "Simulate example detector: CGD" begin
        sim = Simulation{T}(SSD_examples[:CGD])
        simulate!(sim, max_refinements = 2, verbose = true)
        evt = Event([CartesianPoint{T}(5e-3,0,0)])
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-2 )
    end
    @testset "Simulate example detector: Spherical" begin
        sim = Simulation{T}(SSD_examples[:Spherical])
        simulate!(sim, max_refinements = 1, verbose = true)
        evt = Event([CartesianPoint{T}(0,0,0)])
        simulate!(evt, sim)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-4 )
    end 
    @testset "Simulate example detector: Toroidal" begin
        sim = Simulation{T}(SSD_examples[:CoaxialTorus])
        SolidStateDetectors.apply_initial_state!(sim, ElectricPotential)
        simulate!(sim, max_refinements = 1, verbose = true)
        evt = Event([CartesianPoint{T}(0.01,0,0.003)])
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 2e-2 )
    end
    @testset "Simulate example detector: SigGen PPC" begin
        sim = Simulation{T}(SSD_examples[:SigGen])
        simulate!(sim, max_refinements = 1, verbose = true)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 40e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 2e-2 )
    end
end

@testset "ADLChargeDriftModel" begin
    include("ADLChargeDriftModel.jl")
end

# include("ConstructiveSolidGeometry/CSG_test.jl")
include("ConstructiveSolidGeometry/CSG_IO.jl")
# include("ConstructiveSolidGeometry/CSG_decomposition.jl")