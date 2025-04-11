# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test
using SolidStateDetectors

abstract type pseudoGPUArray <: SolidStateDetectors.GPUArrays.AbstractGPUArray{T where T, N where N} end

@testset "Simulate example detector: Inverted Coax" begin
    sim = Simulation{T}(SSD_examples[:InvertedCoax])

    calculate_electric_potential!( 
        sim, 
        device_array_type = pseudoGPUArray, 
        convergence_limit = 1e-5,
        refinement_limits = [0.2, 0.1],
        depletion_handling = true
    )    

    W = SolidStateDetectors.calculate_stored_energy(sim)
    @test isapprox(ustrip(W), 2.41e-4, atol = 1e-5)        
end

@testset "Simulate example detector: CGD" begin
    sim = Simulation{T}(SSD_examples[:CGD])

    calculate_electric_potential!( 
        sim, 
        device_array_type = pseudoGPUArray, 
        convergence_limit = 1e-5,
        refinement_limits = [0.2],
        depletion_handling = true
    )    

    W = SolidStateDetectors.calculate_stored_energy(sim)
    @test isapprox(ustrip(W), 4.39e-6, atol = 1e-5)        
end


