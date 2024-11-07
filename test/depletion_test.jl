using SolidStateDetectors
using Test
using Unitful

T = Float32

@testset "Test depletion estimation" begin
    sim = Simulation{T}("BEGe_01.yaml")
    calculate_electric_potential!(sim, refinement_limits=0.01)
    id = SolidStateDetectors.determine_bias_voltage_contact_id(sim.detector)
    calculate_weighting_potential!(sim, id, refinement_limits=0.01)
    SolidStateDetectors._adapt_weighting_potential_to_electric_potential_grid!(sim, id)
    U_est = ustrip(estimate_depletion_voltage(sim)) # around 2600
    ΔU = 50
    # simulate over and under depletion voltage
    U₋ = U_est - ΔU
    U₊ = U_est + ΔU
    sim.detector = SolidStateDetector(sim.detector, contact_id=1, contact_potential=U₊)
    calculate_electric_potential!(sim, refinement_limits=0.01)
    undepleted = is_depleted(sim.point_types)
    sim.detector = SolidStateDetector(sim.detector, contact_id=1, contact_potential=U₋)
    calculate_electric_potential!(sim, refinement_limits=0.01)
    depleted = is_depleted(sim.point_types)
    @test undepleted && depleted
end
