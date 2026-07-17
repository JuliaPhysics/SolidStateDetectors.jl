# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using SolidStateDetectors
using Test

using SolidStateDetectors: DiscreteAxis, CylindricalGrid, PointTypes, update_bit,
    get_electric_field_from_potential, searchsortednearest

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

@testset "Periodic axis searchsortednearest wraps at the seam" begin
    ax = DiscreteAxis(0.0, 2π, :periodic, :periodic, :closed, :open, collect(range(0.0, 2π, length = 37))[1:end-1])
    @test searchsortednearest(ax, 0.01) == 1
    @test searchsortednearest(ax, 2π - 0.01) == 1       # closest through the wrap
    @test searchsortednearest(ax, 2π + 0.01) == 1
    @test searchsortednearest(ax, -0.01) == 1
    @test searchsortednearest(ax, ax.ticks[end] + 0.01) == 36
end

@testset "Periodic BC on the compressed red-black dimension" begin
    nr, nφ, nz = 3, 4, 6
    axr = DiscreteAxis(0.0, 1.0, :r0, :infinite, :closed, :closed, collect(range(0.0, 1.0, length = nr)))
    axφ = DiscreteAxis(0.0, 2π, :periodic, :periodic, :closed, :open, collect(range(0.0, 2π, length = nφ + 1))[1:end-1])
    axz = DiscreteAxis(0.0, 1.0, :periodic, :periodic, :closed, :open, collect(range(0.0, 1.0, length = nz + 1))[1:end-1])
    grid = CylindricalGrid{Float64}((axr, axφ, axz))
    a = [1000.0 * ir + 100.0 * iφ + iz for ir in 1:nr, iφ in 1:nφ, iz in 1:nz]
    rb = SolidStateDetectors.RBExtBy2Array(a, grid)
    for rbi in (SolidStateDetectors.rb_even, SolidStateDetectors.rb_odd)
        SolidStateDetectors.apply_boundary_conditions_on_x_axis!(rb, rbi, axz, axz.interval, (1.0, 1.0))
    end
    for ir in 1:nr, iφ in 1:nφ
        # the ghost below real z=1 must hold real z=nz of the matching color, and vice versa
        c_zN = iseven(nz + iφ + ir) ? SolidStateDetectors.rb_even : SolidStateDetectors.rb_odd
        c_z1 = iseven(1 + iφ + ir) ? SolidStateDetectors.rb_even : SolidStateDetectors.rb_odd
        @test rb[1, iφ + 1, ir + 1, c_zN] == a[ir, iφ, nz]
        @test rb[end, iφ + 1, ir + 1, c_z1] == a[ir, iφ, 1]
    end
end
