# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).
using Test

using SolidStateDetectors

using Unitful
ϵ0 = SolidStateDetectors.ϵ0 * u"F / m"

T = Float32

@testset "Test real detectors" begin
    @testset "Simulate example detector: Inverted Coax" begin
        sim = Simulation(SSD_examples[:InvertedCoax])
        simulate!(sim, max_refinements = 1, verbose = true)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 40e-3 )]))
        simulate!(evt, sim)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
    @testset "Simulate example detector: Coax" begin
        sim = Simulation(SSD_examples[:Coax])
        simulate!(sim, max_refinements = 0, verbose = true)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(30), 12e-3 )]))
        simulate!(evt, sim)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @test isapprox( signalsum, T(2), atol = 5e-2 )
    end
    @testset "Simulate example detector: BEGe" begin
        sim = Simulation(SSD_examples[:BEGe])
        simulate!(sim, max_refinements = 0, verbose = true)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 20e-3 )]))
        simulate!(evt, sim)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @test isapprox( signalsum, T(2), atol = 1e-2 )
    end
    @testset "Simulate example detector: CGD" begin
        sim = Simulation(SSD_examples[:CGD])
        simulate!(sim, max_refinements = 0, verbose = true)
        evt = Event([CartesianPoint{T}(5e-3, 5e-3, 5e-3)])
        simulate!(evt, sim)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @test isapprox( signalsum, T(2), atol = 1e-2 )
    end
    @testset "Simulate example detector: Spherical" begin
        sim = Simulation(SSD_examples[:Spherical])
        simulate!(sim, max_refinements = 1, verbose = true)
        evt = Event([CartesianPoint{T}(0,0,0)])
        simulate!(evt, sim)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @test isapprox( signalsum, T(2), atol = 1e-2 )
    end 
    @testset "Simulate example detector: SigGen PPC" begin
        sim = Simulation(SSD_examples[:SigGen])
        simulate!(sim, max_refinements = 0, verbose = true)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 40e-3 )]))
        simulate!(evt, sim)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(evt.waveforms[i].value[end])
        end
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
end

@testset "Comparison to analytic solutions" begin
    @testset "Infinite Parallel Plate Capacitor" begin
        sim = Simulation(SSD_examples[:InfiniteParallelPlateCapacitor])
        calculate_electric_potential!(sim, 
            init_grid_spacing = T.( (1e-4, 1e-3, 1e-3) ), 
            max_refinements = 1,
        )
        calculate_electric_field!(sim)
        BV_true = SolidStateDetectors._get_abs_bias_voltage(sim.detector) 
        Δd = (sim.detector.contacts[2].geometry.x[2] - sim.detector.contacts[1].geometry.x[2]) * u"m"
        Δx = (sim.detector.world.intervals[2].right - sim.detector.world.intervals[2].left) * u"m"
        A = Δx * Δx 
        V = A * Δd
        E_true = BV_true / Δd
        ϵr = sim.detector.semiconductors[1].material.ϵ_r
        W_true = uconvert(u"J", (ϵr * ϵ0 * E_true^2 / 2) * A * Δd)
        W_ssd = SolidStateDetectors.calculate_stored_energy(sim);
        C_true = uconvert(u"pF", 2 * W_true / (BV_true^2))
        C_ssd = SolidStateDetectors.calculate_capacitance(sim) 
        @test isapprox(C_ssd, C_true, rtol = 0.001) 
    end

    @testset "InfiniteCoaxialCapacitor" begin
        sim = Simulation(SSD_examples[:InfiniteCoaxialCapacitor])
        calculate_electric_potential!(sim, 
            init_grid_size = (100, 2, 20), 
            max_refinements = 2,
        )
        calculate_electric_field!(sim)
        BV_true = SolidStateDetectors._get_abs_bias_voltage(sim.detector) 
        R1 = sim.detector.contacts[1].geometry.r_interval.right * u"m"
        R2 = sim.detector.contacts[2].geometry.r_interval.left * u"m"
        ϵr = sim.detector.semiconductors[1].material.ϵ_r
        L =  (sim.detector.world.intervals[3].right - sim.detector.world.intervals[3].left) * u"m"
        V = π * R2^2 * L

        C_true = uconvert(u"pF", 2π * ϵr * ϵ0 / log(R2/R1) * L )
        W_true = uconvert(u"J", C_true * BV_true^2 / 2)

        W_ssd = SolidStateDetectors.calculate_stored_energy(sim); 
        C_ssd = SolidStateDetectors.calculate_capacitance(sim);

        @test isapprox(C_ssd, C_true, rtol = 0.01) 
    end
end
