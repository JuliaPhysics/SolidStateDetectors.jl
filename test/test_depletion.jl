using SolidStateDetectors
using Test
using Unitful

T = Float32

@testset "Test depletion estimation" begin
    sim = Simulation{T}(joinpath(@__DIR__, "test_config_files/BEGe_01.yaml"))
    timed_calculate_electric_potential!(sim, refinement_limits=0.01)
    id = SolidStateDetectors.determine_bias_voltage_contact_id(sim.detector)
    timed_calculate_weighting_potential!(sim, id, refinement_limits=0.01)
    SolidStateDetectors._adapt_weighting_potential_to_electric_potential_grid!(
        sim, id)
    U_est = timed_estimate_depletion_voltage(sim) # around 2600
    ΔU = 50u"V"
    # simulate over and under depletion voltage
    U₋ = U_est - ΔU
    U₊ = U_est + ΔU
    sim.detector = SolidStateDetector(sim.detector, contact_id=id, contact_potential=U₊)
    timed_calculate_electric_potential!(sim, refinement_limits=0.01, depletion_handling=true)
    undepleted = !is_depleted(sim.point_types)
    sim.detector = SolidStateDetector(sim.detector, contact_id=id, contact_potential=U₋)
    timed_calculate_electric_potential!(sim, refinement_limits=0.01, depletion_handling=true)
    depleted = is_depleted(sim.point_types)
    @test undepleted && depleted

    # Pass a searching range (with units)
    U_alt = timed_estimate_depletion_voltage(sim, U_est * 1.5, 0u"V", tolerance = 0.1u"V")
    @test abs(U_est - U_alt) < 5u"V"

    @test_throws Exception estimate_depletion_voltage(sim, -abs(U_est), abs(U_est))
    @test_throws Exception estimate_depletion_voltage(sim, -10, 0, tolerance = 20)
    @test_throws Exception estimate_depletion_voltage(sim, 0u"kg", 20u"kg")
    @test_logs (:info,) (:info,) (:warn, r".*not in the specified range.*") estimate_depletion_voltage(sim, U_est/3, 0)
end
