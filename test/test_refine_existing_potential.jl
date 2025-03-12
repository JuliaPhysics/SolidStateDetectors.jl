using SolidStateDetectors
using Test

T = Float32

@testset "Test simulating in two refinement steps" begin
    sim = Simulation{T}("BEGe_01.yaml")
    calculate_electric_potential!(sim, refinement_limits = [0.2, 0.1], depletion_handling = true)

    sim2 = Simulation{T}("BEGe_01.yaml")
    calculate_electric_potential!(sim2, refinement_limits = [0.2], depletion_handling = true)
    calculate_electric_potential!(sim2, refinement_limits = [0.1], depletion_handling = true, initialize = false)

    @test size(sim.electric_potential.data) == size(sim2.electric_potential.data)

    calculate_weighting_potential!(sim, 1, refinement_limits = [0.2, 0.1], depletion_handling = true)
    calculate_weighting_potential!(sim2, 1, refinement_limits = [0.2], depletion_handling = true)
    calculate_weighting_potential!(sim2, 1, refinement_limits = [0.1], depletion_handling = true, initialize = false)

    @test size(sim.weighting_potentials[1].data) == size(sim2.weighting_potentials[1].data)
end
