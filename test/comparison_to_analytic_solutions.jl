ϵ0 = SolidStateDetectors.ϵ0 * u"F / m"
e = SolidStateDetectors.elementary_charge * u"C"

@testset "Infinite Parallel Plate Capacitor" begin
    sim = Simulation{T}(SSD_examples[:InfiniteParallelPlateCapacitor])
    calculate_electric_potential!(sim, 
        init_grid_spacing = T.( (1e-4, 1e-3, 1e-3) ), 
        max_refinements = 1,
    )
    calculate_electric_field!(sim)
    BV_true = SolidStateDetectors._get_abs_bias_voltage(sim.detector) 
    Δd = (sim.detector.contacts[2].geometry.x[2] - sim.detector.contacts[1].geometry.x[2]) * u"m"
    Δx = (sim.detector.world.intervals[2].right - sim.detector.world.intervals[2].left) * u"m"
    A = Δx * Δx 
    V = A * Δd
    E_true = BV_true / Δd
    ϵr = sim.detector.semiconductors[1].material.ϵ_r
    W_true = uconvert(u"J", (ϵr * ϵ0 * E_true^2 / 2) * A * Δd)
    W_ssd = SolidStateDetectors.calculate_stored_energy(sim);
    C_true = uconvert(u"pF", 2 * W_true / (BV_true^2))
    C_ssd = SolidStateDetectors.calculate_capacitance(sim) 
    @testset "Capacity" begin
        @test isapprox(C_ssd, C_true, rtol = 0.001) 
    end
end

struct DummyChargeDensityModel{T} <: SolidStateDetectors.AbstractChargeDensityModel{T} end


@testset "InfiniteCoaxialCapacitor" begin
    sim_cyl = Simulation{T}(SSD_examples[:InfiniteCoaxialCapacitor])
    sim_car = Simulation{T}(SSD_examples[:InfiniteCoaxialCapacitorCartesianCoords])
    BV_true = SolidStateDetectors._get_abs_bias_voltage(sim_cyl.detector)
    ϵr = sim_cyl.detector.semiconductors[1].material.ϵ_r
    R1 = sim_cyl.detector.contacts[1].geometry.r_interval.right * u"m"
    R2 = sim_cyl.detector.contacts[2].geometry.r_interval.left * u"m"
    V1 = sim_cyl.detector.contacts[1].potential * u"V"
    V2 = sim_cyl.detector.contacts[2].potential * u"V"
    L = (sim_cyl.detector.world.intervals[3].right - sim_cyl.detector.world.intervals[3].left) * u"m"
    V = π * R2^2 * L
    intV = (R2^2 - R1^2) * π * L

    calculate_electric_potential!(sim_cyl, 
        init_grid_size = (40, 2, 2), 
        max_refinements = 1,
    )
    calculate_electric_potential!(sim_car, 
        init_grid_size = (40, 40, 2), 
        max_refinements = 1
    )
    calculate_electric_field!(sim_cyl)
    calculate_electric_field!(sim_car)

    C_true = uconvert(u"pF", 2π * ϵr * ϵ0 / log(R2/R1) * L )
    W_true = uconvert(u"J", C_true * BV_true^2 / 2)

    C_cyl_ssd = SolidStateDetectors.calculate_capacitance(sim_cyl)
    C_car_ssd = SolidStateDetectors.calculate_capacitance(sim_car)
    W_cyl_ssd = SolidStateDetectors.calculate_stored_energy(sim_cyl)
    W_car_ssd = SolidStateDetectors.calculate_stored_energy(sim_car)

    @testset "Capacity" begin
        @test isapprox(C_cyl_ssd, C_true, rtol = 0.01) 
        @test isapprox(C_car_ssd, C_true, rtol = 0.06) 
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
    pot_analytic = map(r -> potential_analytic(r), rs_analytic)

    function SolidStateDetectors.get_charge_density(cdm::DummyChargeDensityModel{T}, pt::CylindricalPoint{T})::T where {T <: SSDFloat}
        if ustrip(R1) <= pt[1] <= ustrip(R2)
            return -ustrip(ρ1 / e) * pt[1]^2 
        else 
            return 0
        end
    end

    function SolidStateDetectors.get_charge_density(cdm::DummyChargeDensityModel{T}, pt::CartesianPoint{T})::T where {T <: SSDFloat}
        SolidStateDetectors.get_charge_density(cdm, CylindricalPoint(pt))
    end

    sim_cyl.detector.semiconductors[1].charge_density_model = DummyChargeDensityModel{T}()
    sim_car.detector.semiconductors[1].charge_density_model = DummyChargeDensityModel{T}()

    calculate_electric_potential!(sim_cyl, 
        init_grid_size = (40, 2, 2), 
        max_refinements = 1
    )
    calculate_electric_potential!(sim_car, 
        init_grid_size = (40, 40, 2), 
        max_refinements = 1
    )

    idxR1 = searchsortedfirst( sim_cyl.electric_potential.grid.axes[1], ustrip(R1))
    rs = sim_cyl.electric_potential.grid.axes[1].ticks[idxR1:end]
    potential_rms_cyl = sqrt(sum((ustrip.(map(r -> potential_analytic(r * u"m"), rs)) .- sim_cyl.electric_potential.data[idxR1:end, 1, 1]).^2))
    idxR1 = searchsortedfirst( sim_car.electric_potential.grid.axes[1], ustrip(R1))
    rs = sim_car.electric_potential.grid.axes[1].ticks[idxR1:end]
    idy0 = searchsortedfirst( sim_car.electric_potential.grid.axes[2], 0)
    potential_rms_car = sqrt(sum((ustrip.(map(r -> potential_analytic(r * u"m"), rs)) .- sim_car.electric_potential.data[idxR1:end, idy0, 1]).^2))

    @testset "Potential RMS" begin
        @test potential_rms_cyl < 0.01
        @test potential_rms_car < 1.1
    end
end
