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
    U_est = estimate_depletion_voltage(sim) # around 2600
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
end

@testset "_estimate_depletion_voltage_factor" begin
    # Test functions accept AbstractVector
    center_only = 1.0
    center_zero = -2.0
    neighbors_only = [0.5, 0.0, -0.5]
    neighbors_zero = [-1.0, -2.0, -3.0]
    est = SolidStateDetectors._estimate_depletion_voltage_factor(center_only, center_zero, neighbors_only, neighbors_zero)
    @test est isa Float64
    @test SolidStateDetectors._replace_NaN_with_minimum([NaN, 2.0, 3.0]) == [2.0, 2.0, 3.0]
end
