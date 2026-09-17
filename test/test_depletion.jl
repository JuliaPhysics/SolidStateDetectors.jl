using SolidStateDetectors
using Test
using Unitful

T = Float32

@testset "Test depletion estimation" begin
    sim = Simulation{T}("BEGe_01.yaml")
    calculate_electric_potential!(sim, refinement_limits=0.01)
    id = SolidStateDetectors.determine_bias_voltage_contact_id(sim.detector)
    calculate_weighting_potential!(sim, id, refinement_limits=0.01)
    SolidStateDetectors._adapt_weighting_potential_to_electric_potential_grid!(
        sim, id)
    U_est = estimate_depletion_voltage(sim) # around 2600
    ΔU = 50u"V"
    # simulate over and under depletion voltage
    U₋ = U_est - ΔU
    U₊ = U_est + ΔU
    sim.detector = SolidStateDetector(sim.detector, contact_id=id, contact_potential=U₊)
    calculate_electric_potential!(sim, refinement_limits=0.01, depletion_handling=true)
    undepleted = !is_depleted(sim.point_types)
    sim.detector = SolidStateDetector(sim.detector, contact_id=id, contact_potential=U₋)
    calculate_electric_potential!(sim, refinement_limits=0.01, depletion_handling=true)
    depleted = is_depleted(sim.point_types)
    @test undepleted && depleted

    # Pass a searching range
    U_alt = timed_estimate_depletion_voltage(sim, T(ustrip(u"V", U_est * 1.5)), T(0.0), tolerance = 0.1)
    @test abs(U_est - U_alt) < 5u"V"

    # The analytic local-extremum scan and the bisection fallback should agree on the depletion voltage.
    # Build the same inputs `estimate_depletion_voltage` uses internally.
    simDV = deepcopy(sim)
    SolidStateDetectors._adapt_weighting_potential_to_electric_potential_grid!(simDV, id)
    ϕV = simDV.weighting_potentials[id].data
    ϕρ = simDV.electric_potential.data .- simDV.detector.contacts[id].potential .* ϕV
    pt = simDV.point_types.data
    bulk = findall(pt .& SolidStateDetectors.bulk_bit .> 0)
    inside = findall(pt .& SolidStateDetectors.pn_junction_bit .> 0)
    Umin, Umax = minmax(zero(T), T(1.5 * ustrip(u"V", U_est)))
    U_cand = filter(u -> Umin <= u <= Umax, SolidStateDetectors._find_depletion_voltage_candidates(ϕρ, ϕV, bulk))
    @test length(U_cand) == 1
    U_bis = SolidStateDetectors._find_depletion_voltage_by_bisection(ϕρ, ϕV, inside, bulk, Umin, Umax, T(0.1))
    @test abs(only(U_cand) - U_bis) < 5
    @test abs(only(U_cand) - ustrip(u"V", U_est)) < 5

    @test_throws Exception estimate_depletion_voltage(sim, T(ustrip(u"V", -abs(U_est))), T(ustrip(u"V", abs(U_est))))
    @test_throws Exception estimate_depletion_voltage(sim, T(-10), T(0), tolerance = T(20))
    @test_throws ArgumentError estimate_depletion_voltage(sim, T(ustrip(u"V", U_est/3)), T(0))

    # The depletion voltage is linear in the impurity density: doubling the impurities should double it.
    # Use a fresh simulation so `sim` keeps its state for the tests below. The search range is widened
    # since the doubled depletion voltage lies outside the default range given by the contact potentials.
    sim_2x = Simulation{T}(joinpath(@__DIR__, "BEGe_01.yaml"))
    sim_2x.detector = SolidStateDetector(sim_2x.detector, 2 * sim_2x.detector.semiconductor.impurity_density_model)
    timed_calculate_electric_potential!(sim_2x, refinement_limits=0.01)
    U_est_2x = timed_estimate_depletion_voltage(sim_2x, T(ustrip(u"V", 3 * U_est)), T(0), check_for_depletion = false)
    @test isapprox(U_est_2x, 2 * U_est, rtol = 0.02)
    # The depletion voltage is linear in the impurity density: doubling the impurities should double it.
    # Use a fresh simulation so `sim` keeps its state for the tests below. The search range is widened
    # since the doubled depletion voltage lies outside the default range given by the contact potentials.
    sim_2x = Simulation{T}(joinpath(@__DIR__, "BEGe_01.yaml"))
    sim_2x.detector = SolidStateDetector(sim_2x.detector, 2 * sim_2x.detector.semiconductor.impurity_density_model)
    timed_calculate_electric_potential!(sim_2x, refinement_limits=0.01)
    U_est_2x = timed_estimate_depletion_voltage(sim_2x, T(ustrip(u"V", 3 * U_est)), T(0), check_for_depletion = false)
    @test isapprox(U_est_2x, 2 * U_est, rtol = 0.02)

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
    sim_fresh = Simulation{T}(joinpath(@__DIR__, "BEGe_01.yaml"))
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
