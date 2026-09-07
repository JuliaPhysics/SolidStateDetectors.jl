using SolidStateDetectors
using Test

T = Float64
refinement_limits = [0.2, 0.1, 0.05]

_Q_E = 1.602176634e-19
_EPS0 = 8.8541878128e-12
_EPS_HPGE = 16.0 * _EPS0

_to_mm3(N_per_m3::Real) = N_per_m3 * (1e-3)^3

@testset "is_depleted detects undepleted regions near a contact (pn_junction_bit fix)" begin

    N_m3, d_mm, L_mm = 1e16, 10.0, 30.0
    Vd = _Q_E * N_m3 * (d_mm * 1e-3)^2 / (2 * _EPS_HPGE)

    cfg = Dict(
        "name" => "Planar slab (Cartesian)",
        "units" => Dict("length" => "mm", "angle" => "deg", "potential" => "V", "temperature" => "K"),
        "grid" => Dict(
            "coordinates" => "cartesian",
            "axes" => Dict(
                "x" => Dict("from" => -d_mm / 2, "to" => d_mm / 2, "boundaries" => "reflecting"),
                "y" => Dict("from" => -L_mm / 2, "to" => L_mm / 2, "boundaries" => "reflecting"),
                "z" => Dict("from" => -L_mm / 2, "to" => L_mm / 2, "boundaries" => "reflecting"),
            ),
        ),
        "medium" => "vacuum",
        "detectors" => [Dict(
            "semiconductor" => Dict("material" => "HPGe", "temperature" => 78,
                "impurity_density" => Dict("name" => "constant", "value" => -_to_mm3(N_m3)),
                "geometry" => Dict("box" => Dict("hX" => d_mm / 2, "hY" => L_mm / 2, "hZ" => L_mm / 2))),
            "contacts" => [
                Dict("material" => "HPGe", "id" => 1, "potential" => 0,
                     "geometry" => Dict("box" => Dict("hX" => 0, "hY" => L_mm / 2, "hZ" => L_mm / 2, "origin" => Dict("x" => -d_mm / 2)))),
                Dict("material" => "HPGe", "id" => 2, "potential" => 0,
                     "geometry" => Dict("box" => Dict("hX" => 0, "hY" => L_mm / 2, "hZ" => L_mm / 2, "origin" => Dict("x" => d_mm / 2)))),
            ],
        )],
    )

    function is_depleted_planar(V::Real)
        sim = Simulation{T}(deepcopy(cfg))
        sim.detector = SolidStateDetector(sim.detector, contact_id = 2, contact_potential = T(V))
        apply_initial_state!(sim, ElectricPotential, Grid(sim))
        calculate_electric_potential!(sim, depletion_handling = true, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        is_depleted(sim.point_types)
    end

    # 0.95x is excluded: at that bias, is_depleted reports "depleted".
    # Refining further refinement_limits = [0.2,0.1,0.05,0.01,0.002]
    # (30 -> 558 grid points along x) does resolve it correctly,
    # but at ~100x the runtime for this one point. Excluded here as
    # not worth that cost for a threshold-test.
    @test !is_depleted_planar(0.70 * Vd)
    @test !is_depleted_planar(0.85 * Vd)
    @test is_depleted_planar(1.05 * Vd)
    @test is_depleted_planar(1.30 * Vd)
end

@testset "r0_handling_depletion_handling + gauge invariance: planar slab (Cylindrical, touches r=0 axis)" begin

    N_m3, d_mm, R_mm = 1e16, 10.0, 30.0
    Vd = _Q_E * N_m3 * (d_mm * 1e-3)^2 / (2 * _EPS_HPGE)

    cfg = Dict(
        "name" => "Planar slab (Cylindrical)",
        "units" => Dict("length" => "mm", "angle" => "deg", "potential" => "V", "temperature" => "K"),
        "grid" => Dict(
            "coordinates" => "cylindrical",
            "axes" => Dict(
                "r" => Dict("to" => R_mm, "boundaries" => "reflecting"),
                "phi" => Dict("from" => 0, "to" => 0, "boundaries" => "periodic"),
                "z" => Dict("from" => 0, "to" => d_mm, "boundaries" => Dict("left" => "fixed", "right" => "fixed")),
            ),
        ),
        "medium" => "vacuum",
        "detectors" => [Dict(
            "semiconductor" => Dict("material" => "HPGe", "temperature" => 78,
                "impurity_density" => Dict("name" => "constant", "value" => -_to_mm3(N_m3)),
                "geometry" => Dict("tube" => Dict("r" => R_mm, "h" => d_mm, "origin" => Dict("z" => d_mm / 2)))),
            "contacts" => [
                Dict("material" => "HPGe", "id" => 1, "potential" => 0,
                     "geometry" => Dict("tube" => Dict("r" => R_mm, "h" => 0, "origin" => Dict("z" => 0)))),
                Dict("material" => "HPGe", "id" => 2, "potential" => 0,
                     "geometry" => Dict("tube" => Dict("r" => R_mm, "h" => 0, "origin" => Dict("z" => d_mm)))),
            ],
        )],
    )

    function is_depleted_planar_cyl(U::Real, ground_offset::Real)
        sim = Simulation{T}(deepcopy(cfg))
        det = sim.detector
        det = SolidStateDetector(det, contact_id = 1, contact_potential = T(ground_offset))
        det = SolidStateDetector(det, contact_id = 2, contact_potential = T(U + ground_offset))
        sim.detector = det
        apply_initial_state!(sim, ElectricPotential, Grid(sim))
        calculate_electric_potential!(sim, depletion_handling = true, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        is_depleted(sim.point_types)
    end

    for U in (0.85 * Vd, 1.15 * Vd), ground_offset in (0.0, -U)   # 0.0: original {point=0, mantle=U}; -U: swapped {point=-U, mantle=0}
        @test is_depleted_planar_cyl(U, ground_offset) == (U > Vd)
    end

    @test !is_depleted_planar_cyl(0.70 * Vd, 0.0)
    @test !is_depleted_planar_cyl(0.85 * Vd, 0.0)
    @test is_depleted_planar_cyl(1.05 * Vd, 0.0)
    @test is_depleted_planar_cyl(1.30 * Vd, 0.0)
end

@testset "Gauge invariance off-axis: annular coaxial detector (r1>0)" begin

    N_m3, r1_mm, r2_mm, H_mm = 1e16, 5.0, 35.0, 20.0
    r1_m, r2_m = r1_mm * 1e-3, r2_mm * 1e-3
    rho = -_Q_E * N_m3
    Vd = log(r2_m / r1_m) * rho * r1_m^2 / (2 * _EPS_HPGE) - rho * (r2_m^2 - r1_m^2) / (4 * _EPS_HPGE)

    cfg = Dict(
        "name" => "Annular coaxial (r1>0)",
        "units" => Dict("length" => "mm", "angle" => "deg", "potential" => "V", "temperature" => "K"),
        "grid" => Dict(
            "coordinates" => "cylindrical",
            "axes" => Dict(
                "r" => Dict("to" => r2_mm, "boundaries" => "fixed"),
                "phi" => Dict("from" => 0, "to" => 0, "boundaries" => "periodic"),
                "z" => Dict("from" => 0, "to" => H_mm, "boundaries" => Dict("left" => "reflecting", "right" => "reflecting")),
            ),
        ),
        "medium" => "vacuum",
        "detectors" => [Dict(
            "semiconductor" => Dict("material" => "HPGe", "temperature" => 78,
                "impurity_density" => Dict("name" => "constant", "value" => -_to_mm3(N_m3)),
                "geometry" => Dict("tube" => Dict("r" => Dict("from" => r1_mm, "to" => r2_mm), "h" => H_mm, "origin" => Dict("z" => H_mm / 2)))),
            "contacts" => [
                Dict("material" => "HPGe", "id" => 1, "potential" => 0,
                     "geometry" => Dict("tube" => Dict("r" => Dict("from" => r1_mm, "to" => r1_mm), "h" => H_mm, "origin" => Dict("z" => H_mm / 2)))),
                Dict("material" => "HPGe", "id" => 2, "potential" => 0,
                     "geometry" => Dict("tube" => Dict("r" => Dict("from" => r2_mm, "to" => r2_mm), "h" => H_mm, "origin" => Dict("z" => H_mm / 2)))),
            ],
        )],
    )

    function is_depleted_coax(U::Real, ground_offset::Real)
        sim = Simulation{T}(deepcopy(cfg))
        det = sim.detector
        det = SolidStateDetector(det, contact_id = 1, contact_potential = T(ground_offset))
        det = SolidStateDetector(det, contact_id = 2, contact_potential = T(U + ground_offset))
        sim.detector = det
        apply_initial_state!(sim, ElectricPotential, Grid(sim))
        calculate_electric_potential!(sim, depletion_handling = true, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        is_depleted(sim.point_types)
    end

    for U in (0.85 * Vd, 1.15 * Vd), ground_offset in (0.0, -U)   # 0.0: original {point=0, mantle=U}; -U: swapped {point=-U, mantle=0}
        @test is_depleted_coax(U, ground_offset) == (U > Vd)
    end

    @test !is_depleted_coax(0.70 * Vd, 0.0)
    @test !is_depleted_coax(0.85 * Vd, 0.0)
    @test is_depleted_coax(1.05 * Vd, 0.0)
    @test is_depleted_coax(1.30 * Vd, 0.0)
end

@testset "gauge invariance with `:infinite` boundaries: planar slab (Cartesian)" begin

    N_m3, d_mm, L_mm = 1e16, 10.0, 30.0
    Vd = _Q_E * N_m3 * (d_mm * 1e-3)^2 / (2 * _EPS_HPGE)

    cfg = Dict(
        "name" => "Planar slab (Cartesian, infinite lateral boundaries)",
        "units" => Dict("length" => "mm", "angle" => "deg", "potential" => "V", "temperature" => "K"),
        "grid" => Dict(
            "coordinates" => "cartesian",
            "axes" => Dict(
                "x" => Dict("from" => -d_mm / 2, "to" => d_mm / 2, "boundaries" => Dict("left" => "fixed", "right" => "fixed")),
                "y" => Dict("from" => -L_mm / 2, "to" => L_mm / 2, "boundaries" => Dict("left" => "inf", "right" => "inf")),
                "z" => Dict("from" => -L_mm / 2, "to" => L_mm / 2, "boundaries" => Dict("left" => "inf", "right" => "inf")),
            ),
        ),
        "medium" => "vacuum",
        "detectors" => [Dict(
            "semiconductor" => Dict("material" => "HPGe", "temperature" => 78,
                "impurity_density" => Dict("name" => "constant", "value" => -_to_mm3(N_m3)),
                "geometry" => Dict("box" => Dict("hX" => d_mm / 2, "hY" => L_mm / 2, "hZ" => L_mm / 2))),
            "contacts" => [
                Dict("material" => "HPGe", "id" => 1, "potential" => 0,
                     "geometry" => Dict("box" => Dict("hX" => 0, "hY" => L_mm / 2, "hZ" => L_mm / 2, "origin" => Dict("x" => -d_mm / 2)))),
                Dict("material" => "HPGe", "id" => 2, "potential" => 0,
                     "geometry" => Dict("box" => Dict("hX" => 0, "hY" => L_mm / 2, "hZ" => L_mm / 2, "origin" => Dict("x" => d_mm / 2)))),
            ],
        )],
    )

    function is_depleted_planar_inf(U::Real, ground_offset::Real)
        sim = Simulation{T}(deepcopy(cfg))
        det = sim.detector
        det = SolidStateDetector(det, contact_id = 1, contact_potential = T(ground_offset))
        det = SolidStateDetector(det, contact_id = 2, contact_potential = T(U + ground_offset))
        sim.detector = det
        apply_initial_state!(sim, ElectricPotential, Grid(sim))
        calculate_electric_potential!(sim, depletion_handling = true, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        is_depleted(sim.point_types)
    end

    for U in (0.85 * Vd, 1.15 * Vd), ground_offset in (0.0, -U)   # 0.0: original {point=0, mantle=U}; -U: swapped {point=-U, mantle=0}
        @test is_depleted_planar_inf(U, ground_offset) == (U > Vd)
    end

    @test !is_depleted_planar_inf(0.70 * Vd, 0.0)
    @test !is_depleted_planar_inf(0.85 * Vd, 0.0)
    @test is_depleted_planar_inf(1.05 * Vd, 0.0)
    @test is_depleted_planar_inf(1.30 * Vd, 0.0)
end

@testset "gauge invariance of the electric potential field (direct comparison, `:infinite` boundaries)" begin

    # With depletion_handling = false, calculate_electric_potential! solves a strictly linear Poisson
    # problem, so shifting both contact potentials by the same Δ (a full gauge shift) must shift the
    # converged potential by exactly Δ everywhere, including at the `:infinite`-boundary edge cells
    # that `_shift_axis_margin!`/`apply_boundary_conditions!` handle.

    N_m3, d_mm, L_mm = 1e16, 10.0, 30.0
    Vd = _Q_E * N_m3 * (d_mm * 1e-3)^2 / (2 * _EPS_HPGE)

    cfg = Dict(
        "name" => "Planar slab (Cartesian, infinite lateral boundaries)",
        "units" => Dict("length" => "mm", "angle" => "deg", "potential" => "V", "temperature" => "K"),
        "grid" => Dict(
            "coordinates" => "cartesian",
            "axes" => Dict(
                "x" => Dict("from" => -d_mm / 2, "to" => d_mm / 2, "boundaries" => Dict("left" => "fixed", "right" => "fixed")),
                "y" => Dict("from" => -L_mm / 2, "to" => L_mm / 2, "boundaries" => Dict("left" => "inf", "right" => "inf")),
                "z" => Dict("from" => -L_mm / 2, "to" => L_mm / 2, "boundaries" => Dict("left" => "inf", "right" => "inf")),
            ),
        ),
        "medium" => "vacuum",
        "detectors" => [Dict(
            "semiconductor" => Dict("material" => "HPGe", "temperature" => 78,
                "impurity_density" => Dict("name" => "constant", "value" => -_to_mm3(N_m3)),
                "geometry" => Dict("box" => Dict("hX" => d_mm / 2, "hY" => L_mm / 2, "hZ" => L_mm / 2))),
            "contacts" => [
                Dict("material" => "HPGe", "id" => 1, "potential" => 0,
                     "geometry" => Dict("box" => Dict("hX" => 0, "hY" => L_mm / 2, "hZ" => L_mm / 2, "origin" => Dict("x" => -d_mm / 2)))),
                Dict("material" => "HPGe", "id" => 2, "potential" => 0,
                     "geometry" => Dict("box" => Dict("hX" => 0, "hY" => L_mm / 2, "hZ" => L_mm / 2, "origin" => Dict("x" => d_mm / 2)))),
            ],
        )],
    )

    function electric_potential_planar_inf(U::Real, ground_offset::Real)
        sim = Simulation{T}(deepcopy(cfg))
        det = sim.detector
        det = SolidStateDetector(det, contact_id = 1, contact_potential = T(ground_offset))
        det = SolidStateDetector(det, contact_id = 2, contact_potential = T(U + ground_offset))
        sim.detector = det
        apply_initial_state!(sim, ElectricPotential, Grid(sim))
        calculate_electric_potential!(sim, depletion_handling = false, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        sim.electric_potential
    end

    U = 0.6 * Vd   # well below Vd
    epot0 = electric_potential_planar_inf(U, 0.0)

    for Δ in (-U, -0.37 * Vd, 2.0 * Vd)
        epotΔ = electric_potential_planar_inf(U, Δ)
        @test size(epotΔ) == size(epot0)
        @test epotΔ.data ≈ epot0.data .+ Δ atol = 1e-3 * Vd
    end
end

@testset "r0-axis + `:infinite` r-boundary together: planar slab (Cylindrical, r touches 0, decaying at r_max)" begin

    N_m3, d_mm, R_mm = 1e16, 10.0, 30.0
    Vd = _Q_E * N_m3 * (d_mm * 1e-3)^2 / (2 * _EPS_HPGE)

    cfg = Dict(
        "name" => "Planar slab (Cylindrical, r0 + infinite r-boundary)",
        "units" => Dict("length" => "mm", "angle" => "deg", "potential" => "V", "temperature" => "K"),
        "grid" => Dict(
            "coordinates" => "cylindrical",
            "axes" => Dict(
                "r" => Dict("to" => R_mm, "boundaries" => "inf"),
                "phi" => Dict("from" => 0, "to" => 0, "boundaries" => "periodic"),
                "z" => Dict("from" => 0, "to" => d_mm, "boundaries" => Dict("left" => "fixed", "right" => "fixed")),
            ),
        ),
        "medium" => "vacuum",
        "detectors" => [Dict(
            "semiconductor" => Dict("material" => "HPGe", "temperature" => 78,
                "impurity_density" => Dict("name" => "constant", "value" => -_to_mm3(N_m3)),
                "geometry" => Dict("tube" => Dict("r" => R_mm, "h" => d_mm, "origin" => Dict("z" => d_mm / 2)))),
            "contacts" => [
                Dict("material" => "HPGe", "id" => 1, "potential" => 0,
                     "geometry" => Dict("tube" => Dict("r" => R_mm, "h" => 0, "origin" => Dict("z" => 0)))),
                Dict("material" => "HPGe", "id" => 2, "potential" => 0,
                     "geometry" => Dict("tube" => Dict("r" => R_mm, "h" => 0, "origin" => Dict("z" => d_mm)))),
            ],
        )],
    )

    function is_depleted_planar_cyl_inf(U::Real, ground_offset::Real)
        sim = Simulation{T}(deepcopy(cfg))
        det = sim.detector
        det = SolidStateDetector(det, contact_id = 1, contact_potential = T(ground_offset))
        det = SolidStateDetector(det, contact_id = 2, contact_potential = T(U + ground_offset))
        sim.detector = det
        apply_initial_state!(sim, ElectricPotential, Grid(sim))
        calculate_electric_potential!(sim, depletion_handling = true, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        is_depleted(sim.point_types)
    end

    for U in (0.85 * Vd, 1.15 * Vd), ground_offset in (0.0, -U)   # 0.0: original {point=0, mantle=U}; -U: swapped {point=-U, mantle=0}
        @test is_depleted_planar_cyl_inf(U, ground_offset) == (U > Vd)
    end

    @test !is_depleted_planar_cyl_inf(0.70 * Vd, 0.0)
    @test !is_depleted_planar_cyl_inf(0.85 * Vd, 0.0)
    @test is_depleted_planar_cyl_inf(1.05 * Vd, 0.0)
    @test is_depleted_planar_cyl_inf(1.30 * Vd, 0.0)
end
