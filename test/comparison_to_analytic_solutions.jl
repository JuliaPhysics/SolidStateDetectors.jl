ϵ0 = SolidStateDetectors.ϵ0 * u"F / m"
e = SolidStateDetectors.elementary_charge * u"C"


@testset "Infinite Parallel Plate Capacitor" begin
    sim = Simulation{T}(SSD_examples[:InfiniteParallelPlateCapacitor])
    calculate_electric_potential!(sim, device_array_type = device_array_type, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
    calculate_weighting_potential!(sim, 1, device_array_type = device_array_type, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
    calculate_weighting_potential!(sim, 2, device_array_type = device_array_type, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
    # calculate_electric_field!(sim)
    BV_true = (maximum(broadcast(c -> c.potential, sim.detector.contacts)) - minimum(broadcast(c -> c.potential, sim.detector.contacts)))* u"V"
    # Δd = (sim.detector.contacts[2].decomposed_surfaces[1].loc - sim.detector.contacts[1].decomposed_surfaces[1].loc) * u"m"
    Δd = (sim.detector.contacts[2].geometry.origin[1] - sim.detector.contacts[1].geometry.origin[1]) * u"m"
    Δx = (sim.world.intervals[2].right - sim.world.intervals[2].left) * u"m"
    A = Δx * Δx 
    V = A * Δd
    E_true = BV_true / Δd
    ϵr = sim.detector.semiconductor.material.ϵ_r
    W_true = uconvert(u"J", (ϵr * ϵ0 * E_true^2 / 2) * A * Δd)
    C_true = uconvert(u"pF", 2 * W_true / (BV_true^2))
    C_Analytical = [C_true -C_true; -C_true C_true]
    C_ssd = calculate_capacitance_matrix(sim)
    @testset "Capacity" begin
        @test all(isapprox.(C_ssd, C_Analytical, rtol = 0.001))
    end
end

@testset "Two Spheres Capacitor" begin
    sim = Simulation{T}(SSD_examples[:TwoSpheresCapacitor])
    simulate!(sim, device_array_type = device_array_type, refinement_limits = [0.2, 0.1, 0.05, 0.02], convergence_limit = 1e-7)
    C_ssd = calculate_capacitance_matrix(sim)

    R_1 = sim.detector.contacts[1].geometry.r * u"m"
    R_2 = sim.detector.contacts[2].geometry.r * u"m"
    d = abs((sim.detector.contacts[1].geometry.origin - sim.detector.contacts[2].geometry.origin)[3]) * u"m"

    # Analytical solutions for the elements of the capacitance matrix:
    function c_ii(r1, r2, d; n = 10)
        r1 = ustrip(uconvert(u"m", r1))
        r2 = ustrip(uconvert(u"m", r2))
        d  = ustrip(uconvert(u"m", d))
        c = 0
        F = 4π*ϵ0 * 1u"m"
        u = acosh((d^2 - r1^2 - r2^2) / (2*r1*r2))
        for i in 0:n
            c += inv(r1*sinh(i*u) + r2*sinh((i + 1)*u))
        end
        uconvert(u"pF", F * r1 * r2 * sinh(u) * c)
    end
    function c_ij(r1, r2, d; n = 10)
        r1 = ustrip(uconvert(u"m", r1))
        r2 = ustrip(uconvert(u"m", r2))
        d  = ustrip(uconvert(u"m", d))
        c = 0
        F = 4π*ϵ0 * 1u"m"
        u = acosh((d^2 - r1^2 - r2^2) / (2*r1*r2))
        for i in 1:n
            c += inv(sinh(i*u))
        end
        uconvert(u"pF", -F * r1 * r2 * sinh(u) * c / d)
    end
    A_c11 = c_ii(R_1, R_2, d)
    A_c22 = c_ii(R_2, R_1, d)
    A_c12 = c_ij(R_1, R_2, d)
    C_Analytical = [A_c11 A_c12; A_c12 A_c22]
    @test all(isapprox.(C_Analytical, C_ssd, rtol = 0.03))
end


struct DummyImpurityDensity{T} <: SolidStateDetectors.AbstractImpurityDensity{T} end


@testset "InfiniteCoaxialCapacitor" begin
    sim_cyl = Simulation{T}(SSD_examples[:InfiniteCoaxialCapacitor])
    sim_car = Simulation{T}(SSD_examples[:InfiniteCoaxialCapacitorCartesianCoords])
    BV_true = (maximum(broadcast(c -> c.potential, sim_cyl.detector.contacts)) - minimum(broadcast(c -> c.potential, sim_cyl.detector.contacts)))* u"V"
    ϵr = sim_cyl.detector.semiconductor.material.ϵ_r
    R1 = sim_cyl.detector.contacts[1].geometry.r[1][1] * u"m"
    R2 = sim_cyl.detector.contacts[2].geometry.r[1][1] * u"m"
    V1 = sim_cyl.detector.contacts[1].potential * u"V"
    V2 = sim_cyl.detector.contacts[2].potential * u"V"
    L = (sim_cyl.world.intervals[3].right - sim_cyl.world.intervals[3].left) * u"m"
    V = π * R2^2 * L
    intV = (R2^2 - R1^2) * π * L

    calculate_electric_potential!(sim_cyl, device_array_type = device_array_type, grid = Grid(sim_cyl, max_tick_distance = (0.3u"cm", 15u"°", 1u"m")),
        convergence_limit = 1e-6, refinement_limits = [0.2, 0.1, 0.05, 0.01], use_nthreads = 1, verbose = false
    )
    calculate_electric_potential!(sim_car, device_array_type = device_array_type, grid = Grid(sim_car, max_tick_distance = (0.3u"cm", 0.3u"cm", 1u"m")),
        convergence_limit = 1e-6, refinement_limits = [0.2, 0.1, 0.05, 0.01], use_nthreads = 1, verbose = false
    )
    calculate_weighting_potential!(sim_cyl, 1, device_array_type = device_array_type, grid = Grid(sim_cyl, max_tick_distance = (0.3u"cm", 15u"°", 1u"m")),
        convergence_limit = 1e-6, refinement_limits = [0.2, 0.1, 0.05, 0.01], use_nthreads = 1, verbose = false
    )
    calculate_weighting_potential!(sim_car, 1, device_array_type = device_array_type, grid = Grid(sim_car, max_tick_distance = (0.3u"cm", 0.3u"cm", 1u"m")),
        convergence_limit = 1e-3, refinement_limits = [0.2], use_nthreads = 1, verbose = false
    )
    calculate_weighting_potential!(sim_cyl, 2, device_array_type = device_array_type, grid = Grid(sim_cyl, max_tick_distance = (0.3u"cm", 15u"°", 1u"m")),
        convergence_limit = 1e-6, refinement_limits = [0.2, 0.1, 0.05, 0.01], use_nthreads = 1, verbose = false
    )
    calculate_weighting_potential!(sim_car, 2, device_array_type = device_array_type, grid = Grid(sim_car, max_tick_distance = (0.3u"cm", 0.3u"cm", 1u"m")),
        convergence_limit = 1e-3, refinement_limits = [0.2], use_nthreads = 1, verbose = false
    )

    C_true = uconvert(u"pF", 2π * ϵr * ϵ0 / log(R2/R1) * L )
    C_Analytical = [C_true -C_true; -C_true C_true]

    C_cyl_ssd = calculate_capacitance_matrix(sim_cyl)
    C_car_ssd = calculate_capacitance_matrix(sim_car)

    @testset "Capacity" begin
        @test all(isapprox.(C_cyl_ssd, C_Analytical, rtol = 0.01))
        @test all(isapprox.(C_car_ssd, C_Analytical, rtol = 0.06))
    end

    # Add impurity density and compare resulting potential to analytic solution
    Q = uconvert(u"C", -1e7 * e * u"cm^-3" * intV)
    Δz = L #(sim_cyl.electric_potential.grid.axes[3].ticks[end] - sim_cyl.electric_potential.grid.axes[3].ticks[1])* u"m"
    # Δz = 0.04f0 * u"m"
    ρ1 = 4Q / ((R2^4 - R1^4)*Δz*2π)
    function ρ_r(r) 
        return ρ1 * r^2 
    end
    a = ρ1 / (ϵ0 * ϵr) 
    c1 = uconvert(u"V", ((a / 16) * (R2^4 - R1^4) - (V2-V1)) / log(R2/R1))
    c2 = uconvert(u"V", V2 + c1 * log(ustrip(R2)) - a * R2^4 / 16)  # unit: Volt
    function potential_analytic(r)
        return uconvert(u"V", a * r^4 / 16 - c1 * log(ustrip(r)) + c2)
    end
    rs_analytic = range(R1, stop = R2, length = 500);
    pot_analytic = map(r -> potential_analytic(r), rs_analytic);

    function SolidStateDetectors.get_impurity_density(cdm::DummyImpurityDensity{T}, pt::CylindricalPoint{T})::T where {T}
        if ustrip(R1) <= pt[1] <= ustrip(R2)
            return -ustrip(ρ1 / e) * pt[1]^2 
        else 
            return 0
        end
    end

    function SolidStateDetectors.get_impurity_density(cdm::DummyImpurityDensity{T}, pt::CartesianPoint{T})::T where {T}
        SolidStateDetectors.get_impurity_density(cdm, CylindricalPoint(pt))
    end
    sim_cyl.detector = SolidStateDetector(sim_cyl.detector, DummyImpurityDensity{T}());
    sim_car.detector = SolidStateDetector(sim_car.detector, DummyImpurityDensity{T}());

    calculate_electric_potential!(sim_cyl, device_array_type = device_array_type, grid = Grid(sim_cyl, max_tick_distance = (1u"cm", 15u"°", 1u"m")),
        convergence_limit = 1e-7, refinement_limits = [0.2, 0.1, 0.05], use_nthreads = 1, verbose = false
    )
    calculate_electric_potential!(sim_car, device_array_type = device_array_type, grid = Grid(sim_car, max_tick_distance = (0.5u"cm", 0.5u"cm", 1u"m")),
        convergence_limit = 1e-7, refinement_limits = [0.2, 0.1, 0.05, 0.02], use_nthreads = 1, verbose = false
    )

    idxR1 = searchsortedfirst( sim_cyl.electric_potential.grid.axes[1], ustrip(R1));
    rs = sim_cyl.electric_potential.grid.axes[1].ticks[idxR1:end];
    potential_rms_cyl = sqrt(sum((ustrip.(map(r -> potential_analytic(r * u"m"), rs)) .- sim_cyl.electric_potential.data[idxR1:end, 1, 1]).^2)) / (ustrip(BV_true)*length(rs));
    idxR1 = searchsortedfirst( sim_car.electric_potential.grid.axes[1], ustrip(R1));
    rs = sim_car.electric_potential.grid.axes[1].ticks[idxR1:end];
    idy0 = searchsortedfirst( sim_car.electric_potential.grid.axes[2], 0);
    potential_rms_car = sqrt(sum((ustrip.(map(r -> potential_analytic(r * u"m"), rs)) .- sim_car.electric_potential.data[idxR1:end, idy0, 1]).^2))/ (ustrip(BV_true)*length(rs));

    @testset "Potential RMS" begin
        @test potential_rms_cyl < 1e-3
        @test potential_rms_car < 1e-3
    end
end
