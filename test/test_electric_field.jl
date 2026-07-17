# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using SolidStateDetectors
using Test

using SolidStateDetectors: DiscreteAxis, CylindricalGrid, PointTypes, update_bit,
    get_electric_field_from_potential

@testset "Cylindrical E field at the periodic φ wrap" begin
    T = Float64
    nr, nφ, nz = 6, 36, 5
    rmax, zmax = 0.05, 0.04
    A = 100.0  # V/m
    axr = DiscreteAxis(0.0, rmax, :r0, :infinite, :closed, :closed, collect(range(0.0, rmax, length = nr)))
    axφ = DiscreteAxis(0.0, 2π, :periodic, :periodic, :closed, :open, collect(range(0.0, 2π, length = nφ + 1))[1:end-1])
    axz = DiscreteAxis(0.0, zmax, :infinite, :infinite, :closed, :closed, collect(range(0.0, zmax, length = nz)))
    grid = CylindricalGrid{T}((axr, axφ, axz))

    # V = A * r * cos(φ) = A * x, so E = -∇V = (-A, 0, 0) everywhere
    pot = [A * r * cos(φ) for r in axr.ticks, φ in axφ.ticks, z in axz.ticks]
    epot = ElectricPotential(pot, grid)
    point_types = PointTypes([update_bit for r in axr.ticks, φ in axφ.ticks, z in axz.ticks], grid)

    ef = get_electric_field_from_potential(epot, point_types)
    for iφ in (1, nφ ÷ 2, nφ), ir in (2, nr), iz in (2,)
        E = ef.data[ir, iφ, iz]
        @test E[1] ≈ -A rtol = 0.02
        @test abs(E[2]) < 0.02 * A
        @test abs(E[3]) < 1e-9 * A
    end
end
