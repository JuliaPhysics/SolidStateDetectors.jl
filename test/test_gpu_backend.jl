# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test
using SolidStateDetectors
using Unitful
using JLArrays

T = Float32

# pseudoGPUArray flips `via_KernelAbstractions` so the SOR kernel runs through
# KernelAbstractions, but as an abstract type it survives Adapt unchanged: the
# setup arrays stay `Array`s and none of the `AbstractGPUArray`-dispatched code
# (fused boundary-condition kernels, r0 kernel, device convergence checks) runs.
# The JLArray testsets below cover those paths with a real GPUArray on CPU memory.
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
    @test isapprox(ustrip(W), 4.39e-6, rtol = 0.05)
end

@testset "JLArray backend reproduces the CPU result: Inverted Coax" begin
    simC = Simulation{T}(SSD_examples[:InvertedCoax])
    calculate_electric_potential!(simC, convergence_limit = 1e-5,
        refinement_limits = [0.2, 0.1], depletion_handling = true, verbose = false)
    simJ = Simulation{T}(SSD_examples[:InvertedCoax])
    calculate_electric_potential!(simJ, device_array_type = JLArray, convergence_limit = 1e-5,
        refinement_limits = [0.2, 0.1], depletion_handling = true, verbose = false)
    @test simJ.electric_potential.grid == simC.electric_potential.grid
    @test maximum(abs.(simJ.electric_potential.data .- simC.electric_potential.data)) < T(0.01)
    @test count(simJ.point_types.data .!= simC.point_types.data) <= length(simC.point_types.data) ÷ 10^4
end

@testset "JLArray backend reproduces the CPU result: CGD" begin
    simC = Simulation{T}(SSD_examples[:CGD])
    calculate_electric_potential!(simC, convergence_limit = 1e-5,
        refinement_limits = [0.2], depletion_handling = true, verbose = false)
    simJ = Simulation{T}(SSD_examples[:CGD])
    calculate_electric_potential!(simJ, device_array_type = JLArray, convergence_limit = 1e-5,
        refinement_limits = [0.2], depletion_handling = true, verbose = false)
    @test simJ.electric_potential.grid == simC.electric_potential.grid
    @test maximum(abs.(simJ.electric_potential.data .- simC.electric_potential.data)) < T(0.01)
    @test count(simJ.point_types.data .!= simC.point_types.data) <= length(simC.point_types.data) ÷ 10^4
end


