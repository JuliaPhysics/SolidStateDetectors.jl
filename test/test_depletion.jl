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
    U_est = timed_estimate_depletion_voltage(sim, check_for_depletion = false) # around 2600
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
    @test_throws ArgumentError estimate_depletion_voltage(sim, U_est/3, 0)

    # `adjust_impurity_and_electric_potential_to_match_depletion!` rescales the impurity density and
    # `adjust_bias_and_electric_potential!` swaps in a new contact potential, both
    # analytically (no re-solve) via superposition, so the simulation matches a target
    # depletion voltage `dep` and bias voltage `bias`. When chaining both, match the depletion
    # first and the bias second with `check_against_depletion_voltage = false` (the depletion
    # voltage has just been set, so the check is redundant). Check the round-trip: the estimated
    # depletion voltage should be ≈ `dep` and the bias contact should sit at `bias`.
    dep_target = -2000u"V"   # BEGe_01 depletes at a negative bias (U_est ≈ -2380 V)
    bias_target = -2500u"V"
    imp_model_before = sim.detector.semiconductor.impurity_density_model
    adjust_impurity_and_electric_potential_to_match_depletion!(sim, dep_target, check_for_depletion = false, verbose = false)
    adjust_bias_and_electric_potential!(sim, bias_target, check_against_depletion_voltage = false, verbose = false, reconverge_electric_potential = true)
    @test sim.detector.contacts[id].potential == SolidStateDetectors._parse_value(T, bias_target, SolidStateDetectors.internal_voltage_unit)
    dep_sim = estimate_depletion_voltage(sim, check_for_depletion = false, verbose = false)
    @test abs(dep_sim - dep_target) < 10u"V"
    @test sim.detector.semiconductor.impurity_density_model != imp_model_before

    # Re-run simulation in place and check depletion voltage matches again. This is a check that impurity_density_model and
    # contact_potential where adapted correctly
    timed_calculate_electric_potential!(sim, refinement_limits = 0.01, depletion_handling = true)
    @test abs(estimate_depletion_voltage(sim, check_for_depletion = false, verbose = false) - dep_sim) < 5u"V"

    # Finally, compare to fresh simulation which is changed manually
    sim_fresh = Simulation{T}(joinpath(@__DIR__, "test_config_files/BEGe_01.yaml"))
    sim_fresh.detector = SolidStateDetector(sim_fresh.detector, contact_id = id, contact_potential = bias_target)
    sim_fresh.detector = SolidStateDetector(sim_fresh.detector, sim.detector.semiconductor.impurity_density_model)
    timed_calculate_electric_potential!(sim_fresh, refinement_limits = 0.01, depletion_handling = true)
    @test abs(estimate_depletion_voltage(sim_fresh, check_for_depletion = false, verbose = false) - dep_sim) < 5u"V"

    # Error handling: both functions require the target voltage to share the
    # (non-zero) sign of the relevant reference voltage AND to exceed it in magnitude (the detector
    # must be over-depleted). `sim` now has a depletion voltage ≈ dep_target = -2000 V and a bias of
    # bias_target = -2500 V.

    # `adjust_impurity_and_electric_potential_to_match_depletion!` validates the target depletion voltage
    # against the current bias (-2500 V):
    @test_throws ArgumentError adjust_impurity_and_electric_potential_to_match_depletion!(sim, 2000u"V", verbose = false)    # opposite sign
    @test_throws ArgumentError adjust_impurity_and_electric_potential_to_match_depletion!(sim, -3000u"V", verbose = false)   # |dep| > |bias| (not over-depleted)

    # `adjust_bias_and_electric_potential!` with `check_against_depletion_voltage = true`
    # validates the target bias against the depletion voltage (≈ -2000 V):
    @test_throws ArgumentError adjust_bias_and_electric_potential!(sim, 1000u"V", check_against_depletion_voltage = true, check_for_depletion = false, verbose = false)          # opposite sign
    @test_throws ArgumentError adjust_bias_and_electric_potential!(sim, dep_target / 2, check_against_depletion_voltage = true, check_for_depletion = false, verbose = false)    # |bias| < |dep|
end
