# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test
using SolidStateDetectors
using Unitful

# include("test_utils.jl")

T = Float64
refinement_limits = [0.2, 0.1, 0.05]

_Q_E = 1.602176634e-19
_EPS0 = 8.8541878128e-12
_EPS_HPGE = 16.0 * _EPS0

_to_mm3(N_per_m3::Real) = N_per_m3 * (1e-3)^3

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
        timed_calculate_electric_potential!(sim, depletion_handling = true, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        is_depleted(sim.point_types)
    end

    for U in (0.85 * Vd, 1.02 * Vd, 1.15 * Vd), ground_offset in (0.0, -U)   # 0.0: original {point=0, mantle=U}; -U: swapped {point=-U, mantle=0}
        @test is_depleted_planar_cyl(U, ground_offset) == (U > Vd)
    end

    @test !is_depleted_planar_cyl(0.70 * Vd, 0.0)
    @test !is_depleted_planar_cyl(0.85 * Vd, 0.0)
    @test is_depleted_planar_cyl(1.05 * Vd, 0.0)
    @test is_depleted_planar_cyl(1.30 * Vd, 0.0)

    # Right at the threshold (0.95Vd here), the coarse `refinement_limits` above isn't enough to
    # resolve the shrinking undepleted sliver. A finer grid resolves it correctly.
    function is_depleted_planar_cyl_fine(U::Real, ground_offset::Real)
        sim = Simulation{T}(deepcopy(cfg))
        det = sim.detector
        det = SolidStateDetector(det, contact_id = 1, contact_potential = T(ground_offset))
        det = SolidStateDetector(det, contact_id = 2, contact_potential = T(U + ground_offset))
        sim.detector = det
        apply_initial_state!(sim, ElectricPotential, Grid(sim))
        timed_calculate_electric_potential!(sim, depletion_handling = true, refinement_limits = [refinement_limits; 0.01; 0.002], use_nthreads = 4, verbose = false)
        is_depleted(sim.point_types)
    end
    U = 0.95 * Vd
    for ground_offset in (0.0, -U)
        @test is_depleted_planar_cyl_fine(U, ground_offset) == (U > Vd)
    end

    #  1. depletion_handling = false: a strictly linear Poisson problem, so a full gauge shift Δ
    #     must reproduce epot0 + Δ everywhere. Since `reflecting`/`fixed`/`periodic` have no
    #     absolute reference baked in (unlike `:infinite`), this holds automatically here -- the
    #     point is to confirm the r0-axis handling itself doesn't quietly break that.
    #  2. depletion_handling = true: the same shift claim, with the nonlinear undepleted-region
    #     clamp active. This is the direct probe of `gauge_ref_potential` (the mean-based seed
    #     fix) in this geometry.

    function electric_potential_planar_cyl(U::Real, ground_offset::Real; depletion_handling::Bool)
        sim = Simulation{T}(deepcopy(cfg))
        det = sim.detector
        det = SolidStateDetector(det, contact_id = 1, contact_potential = T(ground_offset))
        det = SolidStateDetector(det, contact_id = 2, contact_potential = T(U + ground_offset))
        sim.detector = det
        apply_initial_state!(sim, ElectricPotential, Grid(sim))
        timed_calculate_electric_potential!(sim, depletion_handling = depletion_handling, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        sim.electric_potential
    end

    U_field = 0.6 * Vd

    for depletion_handling in (false, true)
        epot0 = electric_potential_planar_cyl(U_field, 0.0; depletion_handling)
        for Δ in (-U_field, -0.37 * Vd, 2.0 * Vd)
            epotΔ = electric_potential_planar_cyl(U_field, Δ; depletion_handling)
            @test size(epotΔ) == size(epot0)
            @test epotΔ.data ≈ epot0.data .+ Δ atol = 1e-3 * Vd
        end
    end
end

@testset "Gauge invariance off-axis: annular coaxial detector (r1>0)" begin

    # Two checks here: (1) a simple pass/fail is_depleted(V) check, and (2) a stronger one that
    # predicts where the leftover undepleted spot sits below full depletion, not just whether
    # it's gone. rstar(V) finds that spot by solving Poisson's equation as if the whole volume
    # was already depleted, and locating where that field would be zero. Vd itself is
    # just rstar(V) evaluated where the island shrinks to nothing at the inner contact (r*=r1),
    # so the two formulas are the same physics, not independent checks.

    N_m3, r1_mm, r2_mm, H_mm = 1e16, 5.0, 35.0, 20.0
    r1_m, r2_m = r1_mm * 1e-3, r2_mm * 1e-3
    rho = -_Q_E * N_m3
    Vd = log(r2_m / r1_m) * rho * r1_m^2 / (2 * _EPS_HPGE) - rho * (r2_m^2 - r1_m^2) / (4 * _EPS_HPGE)
    rstar(V) = sqrt((2 * _EPS_HPGE * V / rho + (r2_m^2 - r1_m^2) / 2) / log(r2_m / r1_m))

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
        timed_calculate_electric_potential!(sim, depletion_handling = true, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        is_depleted(sim.point_types)
    end

    for U in (0.85 * Vd, 1.02 * Vd, 1.15 * Vd), ground_offset in (0.0, -U)   # 0.0: original {point=0, mantle=U}; -U: swapped {point=-U, mantle=0}
        @test is_depleted_coax(U, ground_offset) == (U > Vd)
    end

    @test !is_depleted_coax(0.70 * Vd, 0.0)
    @test !is_depleted_coax(0.85 * Vd, 0.0)
    @test is_depleted_coax(1.05 * Vd, 0.0)
    @test is_depleted_coax(1.30 * Vd, 0.0)

    function undepleted_r_range(U::Real, ground_offset::Real)
        sim = Simulation{T}(deepcopy(cfg))
        det = sim.detector
        det = SolidStateDetector(det, contact_id = 1, contact_potential = T(ground_offset))
        det = SolidStateDetector(det, contact_id = 2, contact_potential = T(U + ground_offset))
        sim.detector = det
        apply_initial_state!(sim, ElectricPotential, Grid(sim))
        timed_calculate_electric_potential!(sim, depletion_handling = true, refinement_limits = [refinement_limits; 0.02], use_nthreads = 4, verbose = false)
        pt = sim.point_types.data
        r = sim.point_types.grid.r
        undep = SolidStateDetectors.undepleted_bit
        pnb = SolidStateDetectors.pn_junction_bit
        inact = SolidStateDetectors.inactive_layer_bit
        rs = Float64[]
        for i1 in eachindex(r), i2 in axes(pt, 2), i3 in axes(pt, 3)
            b = pt[i1, i2, i3]
            if (undep & b) > 0 && (pnb & b) > 0 && (inact & b) == 0
                push!(rs, r[i1])
            end
        end
        isempty(rs) ? (NaN, NaN) : (minimum(rs), maximum(rs))
    end

    # Generous margin (3mm) relative to the grid tick spacing near the pinch-off point at this refinement
    margin = 0.003
    for frac in (0.5, 0.7, 0.9), ground_offset_frac in (0.0, -1.0)   # 0.0: original; -1.0: swapped ({point=-U, mantle=0})
        U = frac * Vd
        ground_offset = ground_offset_frac * U
        r_pred = rstar(U)
        rmin, rmax = undepleted_r_range(U, ground_offset)
        @test !isnan(rmin)
        @test rmin - margin <= r_pred <= rmax + margin
    end
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
        timed_calculate_electric_potential!(sim, depletion_handling = true, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
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

    # Two related but distinct claims checked together here, since both reuse the same geometry:
    #  1. depletion_handling = false: a strictly linear Poisson problem, so shifting all
    #     contact potentials by the same Δ must shift the converged potential by exactly Δ
    #     everywhere, including at the `:infinite`-boundary edge cells that
    #     `_shift_axis_margin!`/`apply_boundary_conditions!` handle. Isolates the
    #     `:infinite`-specific Δ/gauge_ref_potential fix.
    #  2. depletion_handling = true: the same shift claim, but with the nonlinear
    #     undepleted-region clamp active. Shifting all contacts by Δ is a pure gauge
    #     transformation -- it can't change the field or the depletion state, only the
    #     absolute reference -- so epotΔ ≈ epot0 + Δ must still hold. If the interior seed
    #     (gauge_ref_potential, the mean-based fix) weren't gauge-consistent, the bistable
    #     clamp could converge to a different undepleted-region shape between the two solves,
    #     showing up here as a mismatch localized near that region, not just a coarse
    #     `is_depleted` disagreement.

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

    function electric_potential_planar_inf(U::Real, ground_offset::Real; depletion_handling::Bool)
        sim = Simulation{T}(deepcopy(cfg))
        det = sim.detector
        det = SolidStateDetector(det, contact_id = 1, contact_potential = T(ground_offset))
        det = SolidStateDetector(det, contact_id = 2, contact_potential = T(U + ground_offset))
        sim.detector = det
        apply_initial_state!(sim, ElectricPotential, Grid(sim))
        timed_calculate_electric_potential!(sim, depletion_handling = depletion_handling, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        sim.electric_potential
    end

    U = 0.6 * Vd   # well below Vd: keeps the dh=false solve well posed and the dh=true solve clear of near-threshold stagnation

    for depletion_handling in (false, true)
        epot0 = electric_potential_planar_inf(U, 0.0; depletion_handling)
        for Δ in (-U, -0.37 * Vd, 2.0 * Vd)
            epotΔ = electric_potential_planar_inf(U, Δ; depletion_handling)
            @test size(epotΔ) == size(epot0)
            @test epotΔ.data ≈ epot0.data .+ Δ atol = 1e-3 * Vd
        end
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
        timed_calculate_electric_potential!(sim, depletion_handling = true, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
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

@testset "gauge invariance under p<->n sign flip: is_depleted threshold (planar slab, Cartesian, reflecting boundaries)" begin

    N_m3, d_mm, L_mm = 1e16, 10.0, 30.0
    Vd = _Q_E * N_m3 * (d_mm * 1e-3)^2 / (2 * _EPS_HPGE)

    function is_depleted_planar_signed(N_sign::Real, V::Real)
        cfg = Dict(
            "name" => "Planar slab (Cartesian, p<->n sign flip)",
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
                    "impurity_density" => Dict("name" => "constant", "value" => N_sign * _to_mm3(N_m3)),
                    "geometry" => Dict("box" => Dict("hX" => d_mm / 2, "hY" => L_mm / 2, "hZ" => L_mm / 2))),
                "contacts" => [
                    Dict("material" => "HPGe", "id" => 1, "potential" => 0,
                         "geometry" => Dict("box" => Dict("hX" => 0, "hY" => L_mm / 2, "hZ" => L_mm / 2, "origin" => Dict("x" => -d_mm / 2)))),
                    Dict("material" => "HPGe", "id" => 2, "potential" => V,
                         "geometry" => Dict("box" => Dict("hX" => 0, "hY" => L_mm / 2, "hZ" => L_mm / 2, "origin" => Dict("x" => d_mm / 2)))),
                ],
            )],
        )
        sim = Simulation{T}(cfg)
        apply_initial_state!(sim, ElectricPotential, Grid(sim))
        timed_calculate_electric_potential!(sim, depletion_handling = true, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        is_depleted(sim.point_types)
    end

    # p-type (N < 0) with a positive bias vs. its exact n-type mirror (N > 0) with the
    # negative bias of the same magnitude: both must agree on depletion at every fraction
    # of |Vd|, not just for the sign of V that happens to match the original convention.
    for frac in (0.70, 0.85, 1.05, 1.30)
        V = frac * Vd
        @test is_depleted_planar_signed(-1, V) == is_depleted_planar_signed(1, -V)
    end

    @test !is_depleted_planar_signed(-1, 0.70 * Vd)
    @test !is_depleted_planar_signed(-1, 0.85 * Vd)
    @test is_depleted_planar_signed(-1, 1.05 * Vd)
    @test is_depleted_planar_signed(-1, 1.30 * Vd)

    @test !is_depleted_planar_signed(1, -0.70 * Vd)
    @test !is_depleted_planar_signed(1, -0.85 * Vd)
    @test is_depleted_planar_signed(1, -1.05 * Vd)
    @test is_depleted_planar_signed(1, -1.30 * Vd)
end

@testset "gauge invariance under p<->n sign flip (direct field comparison, `:infinite` boundaries)" begin

    # With depletion_handling = false, calculate_electric_potential! solves a strictly linear
    # Poisson problem, so negating every contact potential together with the impurity density
    # must negate the converged potential everywhere exactly -- PROVIDED the far boundary is
    # actually far. `:infinite` only ever decays each edge cell toward a single scalar
    # (`minimum_applied_potential`, restored here as the decay target: unlike `mean`,
    # it doesn't perturb the near-boundary field for realistic asymmetric detectors, e.g. the
    # isochrone-test), and that scalar is NOT antisymmetric under this
    # sign flip. This makes `:infinite`'s own truncation edge NOT exactly sign-flip invariant --
    # the truncation edge itself is never going to be exactly sign-flip invariant, at any domain
    # size, that's a property of the boundary formula, not something that decays away with distance.
    #
    # But this is a boundary-truncation artifact, not a bulk-physics one: with L_mm = 30 (only
    # 15mm from contact to truncation) the mismatch measured ~11V (~3% of U) right at the edge
    # and ~0.1V in the interior; enlarging L_mm to 120 (so `:infinite` operates in the regime it
    # actually approximates -- far from the region of interest -- rather than right on top of
    # it) shrinks the near-center mismatch to ~3e-8V, not the edge mismatch itself, which stays
    # nonzero (measured ~0.5V right at the edge cell even at L_mm=120). That's the real fix: give
    # `:infinite` the separation it needs so this residual decays away before reaching the region
    # of interest, not force its reference to something that corrupts other detectors.
    N_m3, d_mm, L_mm = 1e16, 10.0, 120.0
    Vd = _Q_E * N_m3 * (d_mm * 1e-3)^2 / (2 * _EPS_HPGE)

    function electric_potential_planar_inf_signed(N_sign::Real, V::Real)
        cfg = Dict(
            "name" => "Planar slab (Cartesian, infinite lateral boundaries, p<->n sign flip)",
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
                    "impurity_density" => Dict("name" => "constant", "value" => N_sign * _to_mm3(N_m3)),
                    "geometry" => Dict("box" => Dict("hX" => d_mm / 2, "hY" => L_mm / 2, "hZ" => L_mm / 2))),
                "contacts" => [
                    Dict("material" => "HPGe", "id" => 1, "potential" => 0,
                         "geometry" => Dict("box" => Dict("hX" => 0, "hY" => L_mm / 2, "hZ" => L_mm / 2, "origin" => Dict("x" => -d_mm / 2)))),
                    Dict("material" => "HPGe", "id" => 2, "potential" => V,
                         "geometry" => Dict("box" => Dict("hX" => 0, "hY" => L_mm / 2, "hZ" => L_mm / 2, "origin" => Dict("x" => d_mm / 2)))),
                ],
            )],
        )
        sim = Simulation{T}(cfg)
        apply_initial_state!(sim, ElectricPotential, Grid(sim))
        timed_calculate_electric_potential!(sim, depletion_handling = false, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        sim.electric_potential
    end

    U = 0.6 * Vd   # well below Vd, so the (linear, depletion_handling = false) solve is well posed
    epot_p = electric_potential_planar_inf_signed(-1, U)     # p-type, positive bias
    epot_n = electric_potential_planar_inf_signed(1, -U)     # n-type mirror, negative bias

    @test size(epot_n) == size(epot_p)

    y, z = epot_p.grid[2], epot_p.grid[3]
    y_lo, y_hi = findfirst(v -> v > -0.35 * L_mm * 1e-3, y), findlast(v -> v < 0.35 * L_mm * 1e-3, y)
    z_lo, z_hi = findfirst(v -> v > -0.35 * L_mm * 1e-3, z), findlast(v -> v < 0.35 * L_mm * 1e-3, z)
    @test epot_n.data[:, y_lo:y_hi, z_lo:z_hi] ≈ -epot_p.data[:, y_lo:y_hi, z_lo:z_hi] atol = 1e-3 * Vd
end

@testset "gauge invariance under p<->n sign flip: estimate_depletion_voltage matches in magnitude" begin

    # L_mm = 120 (not 30) for the same reason as the direct-field-comparison testset above:
    # `:infinite` only ever decays toward one scalar, and `minimum_applied_potential` (restored
    # as that target -- `mean` corrupts realistic asymmetric detectors, see the isochrone-test
    # regression) isn't antisymmetric under this sign flip, so the truncation edge itself isn't
    # exactly sign-flip invariant at any domain size. But it's a truncation artifact that decays
    # fast with distance from the region of interest so giving `:infinite` the separation
    # it actually needs fixes this, rather than forcing its reference to something that's wrong
    # for other detectors.
    N_m3, d_mm, L_mm = 1e16, 10.0, 120.0
    Vd = _Q_E * N_m3 * (d_mm * 1e-3)^2 / (2 * _EPS_HPGE)

    function dep_voltage_planar_inf_signed(N_sign::Real, V_seed::Real)
        cfg = Dict(
            "name" => "Planar slab (Cartesian, infinite lateral boundaries, p<->n sign flip)",
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
                    "impurity_density" => Dict("name" => "constant", "value" => N_sign * _to_mm3(N_m3)),
                    "geometry" => Dict("box" => Dict("hX" => d_mm / 2, "hY" => L_mm / 2, "hZ" => L_mm / 2))),
                "contacts" => [
                    Dict("material" => "HPGe", "id" => 1, "potential" => 0,
                         "geometry" => Dict("box" => Dict("hX" => 0, "hY" => L_mm / 2, "hZ" => L_mm / 2, "origin" => Dict("x" => -d_mm / 2)))),
                    Dict("material" => "HPGe", "id" => 2, "potential" => V_seed,
                         "geometry" => Dict("box" => Dict("hX" => 0, "hY" => L_mm / 2, "hZ" => L_mm / 2, "origin" => Dict("x" => d_mm / 2)))),
                ],
            )],
        )
        sim = Simulation{T}(cfg)
        apply_initial_state!(sim, ElectricPotential, Grid(sim))
        timed_calculate_electric_potential!(sim, depletion_handling = true, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        id = SolidStateDetectors.determine_bias_voltage_contact_id(sim.detector)
        timed_calculate_weighting_potential!(sim, id, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        SolidStateDetectors._adapt_weighting_potential_to_electric_potential_grid!(sim, id)
        timed_estimate_depletion_voltage(sim, check_for_depletion = false, verbose = false)
    end

    dep_p = dep_voltage_planar_inf_signed(-1, 1.3 * Vd)    # p-type, positive bias
    dep_n = dep_voltage_planar_inf_signed(1, -1.3 * Vd)    # n-type mirror, negative bias

    @test isapprox(ustrip(dep_p), Vd, atol = 0.1 * Vd)
    @test isapprox(ustrip(dep_n), -Vd, atol = 0.1 * Vd)

    @test isapprox(abs(ustrip(dep_p)), abs(ustrip(dep_n)), atol = 1.0)
end

@testset "gauge invariance inside a grounded cryostat (Faraday cage), and why `:infinite` can't be made exactly gauge invariant" begin

    function cryostat_cfg(R_shell::Real, Z_shell::Real; V_point::Real = 0.0, V_mantle::Real = 4000.0)
        point_contact = Dict(
            "name" => "p contact", "id" => 1, "material" => "HPGe", "potential" => V_point,
            "geometry" => Dict("tube" => Dict("r" => 1.5, "h" => 0.1, "translate" => Dict("z" => 0.05))),
        )
        mantle_geom_parts = [
            Dict("tube" => Dict("r" => Dict("from" => 30.0, "to" => 30.0), "h" => 30.0, "origin" => Dict("z" => 15.0))),
            Dict("tube" => Dict("r" => 30.0, "h" => 0.0, "origin" => Dict("z" => 30.0))),
            Dict("tube" => Dict("r" => Dict("from" => 1.5, "to" => 30.0), "h" => 0.0, "origin" => Dict("z" => 0.0))),
        ]
        mantle_contact = Dict(
            "name" => "n contact", "id" => 2, "material" => "HPGe", "potential" => V_mantle,
            "geometry" => Dict("union" => mantle_geom_parts),
        )
        semiconductor = Dict(
            "material" => "HPGe", "bulk_type" => "p", "temperature" => 77.0,
            "impurity_density" => Dict("name" => "constant", "value" => "1e7cm^-3"),
            "geometry" => Dict("tube" => Dict("r" => 30.0, "h" => 30.0, "translate" => Dict("z" => 15.0))),
        )
        shell_parts = [
            Dict("tube" => Dict("r" => Dict("from" => R_shell, "to" => R_shell + 2.0), "h" => 2 * Z_shell, "origin" => Dict("z" => 0.0))),
            Dict("tube" => Dict("r" => Dict("from" => 0.0, "to" => R_shell + 2.0), "h" => 2.0, "origin" => Dict("z" => Z_shell))),
            Dict("tube" => Dict("r" => Dict("from" => 0.0, "to" => R_shell + 2.0), "h" => 2.0, "origin" => Dict("z" => -Z_shell))),
        ]
        cryostat = Dict("name" => "Cryostat", "material" => "Al", "potential" => 0.0, "geometry" => Dict("union" => shell_parts))
        detector = Dict("semiconductor" => semiconductor, "contacts" => [point_contact, mantle_contact], "passives" => [cryostat])
        Dict(
            "name" => "coax enclosed in a grounded cryostat",
            "units" => Dict("length" => "mm", "angle" => "deg", "potential" => "V", "temperature" => "K"),
            "grid" => Dict(
                "coordinates" => "cylindrical",
                "axes" => Dict(
                    "r" => Dict("from" => 0.0, "to" => R_shell + 5.0, "boundaries" => "fixed"),
                    "phi" => Dict("from" => 0, "to" => 0, "boundaries" => "periodic"),
                    "z" => Dict("from" => -Z_shell - 5.0, "to" => Z_shell + 5.0, "boundaries" => "fixed"),
                ),
            ),
            "medium" => "vacuum",
            "detectors" => [detector],
        )
    end

    function cryostat_potential(R_shell::Real, Z_shell::Real; V_point::Real = 0.0, V_mantle::Real = 4000.0)
        sim = Simulation{T}(cryostat_cfg(R_shell, Z_shell; V_point = V_point, V_mantle = V_mantle))
        timed_calculate_electric_potential!(sim, refinement_limits = [0.2, 0.1], use_nthreads = 4, verbose = false, max_n_iterations = 30000)
        sim.electric_potential
    end

    function at(pot, r, z)
        ir = findmin(abs.(pot.grid[1] .- r))[2]
        iz = findmin(abs.(pot.grid[3] .- z))[2]
        pot.data[ir, 1, iz]
    end

    # (1) Shift invariance emerges near the detector as the grounded shell grows.
    R_big, Z_big = 300.0, 300.0
    pot_orig  = cryostat_potential(R_big, Z_big; V_point = 0.0, V_mantle = 4000.0)
    pot_shift = cryostat_potential(R_big, Z_big; V_point = -4000.0, V_mantle = 0.0)

    # near-detector interior point, well inside any of the shell sizes tried below
    r_mid, z_mid = 0.010, 0.010
    v_orig  = at(pot_orig, r_mid, z_mid)
    v_shift = at(pot_shift, r_mid, z_mid)
    @test isapprox(v_shift + 4000, v_orig, atol = 1.0)   # <1V out of a 4000V shift, near the detector

    # (2) ...but NOT right at the shell itself, which is pinned to an absolute reference and
    # cannot float with the contacts -- confirming shift invariance is an interior/emergent
    # property of a large enough grounded enclosure, not a property `:infinite` could ever
    # bake into its truncation edge exactly, since there is no explicit shell there to anchor to.
    r_edge, z_edge = R_big * 1e-3 * 0.9, Z_big * 1e-3 * 0.9
    v_orig_edge  = at(pot_orig, r_edge, z_edge)
    v_shift_edge = at(pot_shift, r_edge, z_edge)
    @test abs((v_shift_edge + 4000) - v_orig_edge) > 100.0   # large, real mismatch near the shell

    # (3) Sample the vacuum right next to the crystal surface -- exactly where `:infinite` would
    # have to guess a reference -- as the distant shell moves farther out (50mm -> 100mm ->
    # 300mm). T
    #  (a) the value is converging to a real limit, not drifting forever -- each step's change is
    #      smaller than the last.
    #  (b) that limit is neither `gauge_ref_potential` (mean of the two contact potentials) nor
    #      `minimum_applied_potential` (their min) -- so no simple statistic of the contact
    #      voltages alone gives the right answer.
    r_gap, z_gap = 0.0325, 0.0325
    v_gap_small  = at(cryostat_potential(50.0, 50.0), r_gap, z_gap)
    v_gap_medium = at(cryostat_potential(100.0, 100.0), r_gap, z_gap)
    v_gap_large  = at(pot_orig, r_gap, z_gap)
    @test abs(v_gap_large - v_gap_medium) < abs(v_gap_medium - v_gap_small)

    mean_candidate = (0.0 + 4000.0) / 2
    min_candidate  = min(0.0, 4000.0)
    @test !isapprox(v_gap_large, mean_candidate, atol = 500.0)
    @test !isapprox(v_gap_large, min_candidate, atol = 500.0)
end

@testset "gauge invariance: is_depleted on the packaged InvertedCoax example detector" begin

    function is_depleted_ivc(U::Real, ground_offset::Real)
        sim = Simulation{T}(SSD_examples[:InvertedCoax])
        det = sim.detector
        det = SolidStateDetector(det, contact_id = 1, contact_potential = T(ground_offset))
        det = SolidStateDetector(det, contact_id = 2, contact_potential = T(U + ground_offset))
        sim.detector = det
        apply_initial_state!(sim, ElectricPotential, Grid(sim))
        timed_calculate_electric_potential!(sim, depletion_handling = true, refinement_limits = refinement_limits, use_nthreads = 4, verbose = false)
        is_depleted(sim.point_types)
    end

    for U in (2012.0, 2014.0, 2016.0, 2018.0, 2020.0, 2022.0)
        @test is_depleted_ivc(U, 0.0) == is_depleted_ivc(U, -U)   # 0.0: {point=0, mantle=U}; -U: {point=-U, mantle=0}
    end

    @test !is_depleted_ivc(2012.0, 0.0)
    @test !is_depleted_ivc(2014.0, 0.0)
    @test is_depleted_ivc(2016.0, 0.0)
    @test is_depleted_ivc(2018.0, 0.0)
    @test is_depleted_ivc(2020.0, 0.0)
    @test is_depleted_ivc(2022.0, 0.0)
end
