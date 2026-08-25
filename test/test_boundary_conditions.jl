# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test
using SolidStateDetectors

import SolidStateDetectors: DiscreteAxis, Interval, KernelAbstractions, rbidx,
    apply_boundary_conditions_on_x_axis!, apply_boundary_conditions_on_y_axis!,
    apply_boundary_conditions_on_z_axis!, apply_boundary_conditions_on_r_axis!,
    apply_boundary_conditions_on_φ_axis!, apply_boundary_conditions_at_r0!,
    r0_bc_kernel!, _ghost_ops, _ghost_ops_compressed

T = Float64

# The red-black potential arrays carry one ghost row per face; the methods under
# test rewrite ghost rows from interior rows of the same parity slab. Distinct
# entries make every source row uniquely identifiable.
fresh_rbpot() = reshape(T.(1:(6 * 7 * 8 * 2)), 6, 7, 8, 2)

cc_ax(bl, br) = DiscreteAxis{T, bl, br}(Interval{:closed, :closed, T}(0, 1), T[0, 1])
co_ax(bl, br) = DiscreteAxis{T, bl, br}(Interval{:closed, :open, T}(0, 1), T[0, 1])

const GBF = (T(0.25), T(0.5))

# expected ghost writes as (side, source_row, factor)
lo(src, fac = one(T)) = (:lo, src, T(fac))
hi(src, fac = one(T)) = (:hi, src, T(fac))

function expected_result(A::Array{T, 4}, dim::Int, writes)
    B = copy(A)
    slab = view(B, :, :, :, 1)
    for (side, src, fac) in writes
        dst = side === :lo ? 1 : size(B, dim)
        selectdim(slab, dim, dst) .= fac .* selectdim(slab, dim, src)
    end
    B
end

function apply_op_pair!(A::Array{T, 4}, dim::Int, ops)
    slab = view(A, :, :, :, 1)
    for (i, (active, src, fac)) in enumerate(ops)
        active || continue
        dst = i == 1 ? 1 : size(A, dim)
        selectdim(slab, dim, dst) .= fac .* selectdim(slab, dim, src)
    end
    A
end

# Ghost-row semantics per (left, right) boundary type. The compressed
# (red-black) dimension has its interior neighbors at rows 2 / end-1; the
# non-compressed dimensions mirror at 3 / end-2 for reflecting faces.
compressed_writes(n) = Dict(
    (:fixed, :fixed)           => (),
    (:infinite, :infinite)     => (lo(2, GBF[1]), hi(n - 1, GBF[2])),
    (:infinite, :fixed)        => (lo(2, GBF[1]),),
    (:fixed, :infinite)        => (hi(n - 1, GBF[2]),),
    (:infinite, :reflecting)   => (lo(2, GBF[1]), hi(n - 1)),
    (:reflecting, :infinite)   => (lo(2), hi(n - 1, GBF[2])),
    (:reflecting, :reflecting) => (lo(2), hi(n - 1)),
    (:fixed, :reflecting)      => (hi(n - 1),),
    (:reflecting, :fixed)      => (lo(2),),
    (:periodic, :periodic)     => (lo(n - 1), hi(2)),
)
noncompressed_writes(n) = Dict(
    (:fixed, :fixed)           => (),
    (:infinite, :infinite)     => (lo(2, GBF[1]), hi(n - 1, GBF[2])),
    (:infinite, :fixed)        => (lo(2, GBF[1]),),
    (:fixed, :infinite)        => (hi(n - 1, GBF[2]),),
    (:infinite, :reflecting)   => (lo(2, GBF[1]), hi(n - 2)),
    (:reflecting, :infinite)   => (lo(3), hi(n - 1, GBF[2])),
    (:reflecting, :reflecting) => (lo(3), hi(n - 2)),
    (:fixed, :reflecting)      => (hi(n - 2),),
    (:reflecting, :fixed)      => (lo(3),),
    (:periodic, :periodic)     => (lo(n - 1), hi(2)),
)
# r-axis outer face reads end-2 even for the infinite decay (whose factor is
# defined from the end-1 tick); pins current behavior, see issue #617
r_axis_writes(n) = Dict(
    (:r0, :infinite)   => (hi(n - 2, GBF[2]),),
    (:r0, :reflecting) => (hi(n - 2),),
    (:r0, :fixed)      => (),
)

axis_for(bl, br) = bl === :periodic ? co_ax(bl, br) : cc_ax(bl, br)

@testset "CPU ghost-row methods per axis and boundary types" begin
    @testset "compressed dimension (Cartesian x / cylindrical z)" begin
        for (key, writes) in compressed_writes(6)
            A = fresh_rbpot()
            ax = axis_for(key...)
            apply_boundary_conditions_on_x_axis!(A, 1, ax, ax.interval, GBF)
            @test A == expected_result(fresh_rbpot(), 1, writes)
        end
    end
    @testset "non-compressed Cartesian y" begin
        for (key, writes) in noncompressed_writes(7)
            A = fresh_rbpot()
            ax = axis_for(key...)
            apply_boundary_conditions_on_y_axis!(A, 1, ax, ax.interval, GBF)
            @test A == expected_result(fresh_rbpot(), 2, writes)
        end
    end
    @testset "non-compressed Cartesian z" begin
        # the mixed infinite/reflecting combinations have no z-axis methods
        # (see issue #617), so only the implemented ones are pinned here
        for (key, writes) in noncompressed_writes(8)
            key in ((:infinite, :reflecting), (:reflecting, :infinite)) && continue
            A = fresh_rbpot()
            ax = axis_for(key...)
            apply_boundary_conditions_on_z_axis!(A, 1, ax, ax.interval, GBF)
            @test A == expected_result(fresh_rbpot(), 3, writes)
        end
    end
    @testset "cylindrical r" begin
        for (key, writes) in r_axis_writes(8)
            A = fresh_rbpot()
            ax = cc_ax(key...)
            apply_boundary_conditions_on_r_axis!(A, 1, ax, ax.interval, GBF)
            @test A == expected_result(fresh_rbpot(), 3, writes)
        end
    end
    @testset "cylindrical φ" begin
        for (key, writes) in ((:periodic, :periodic) => (lo(6), hi(2)),
                              (:reflecting, :reflecting) => (lo(3), hi(5)))
            A = fresh_rbpot()
            ax = axis_for(key...)
            apply_boundary_conditions_on_φ_axis!(A, 1, ax, ax.interval)
            @test A == expected_result(fresh_rbpot(), 2, writes)
        end
    end
end

@testset "GPU ghost-face op tables match the CPU methods" begin
    for (key, _) in compressed_writes(6)
        A = fresh_rbpot()
        ax = axis_for(key...)
        apply_boundary_conditions_on_x_axis!(A, 1, ax, ax.interval, GBF)
        B = apply_op_pair!(fresh_rbpot(), 1, _ghost_ops_compressed(ax, 6, GBF))
        @test A == B
    end
    for (key, _) in noncompressed_writes(7)
        A = fresh_rbpot()
        ax = axis_for(key...)
        apply_boundary_conditions_on_y_axis!(A, 1, ax, ax.interval, GBF)
        B = apply_op_pair!(fresh_rbpot(), 2, _ghost_ops(ax, 7, GBF))
        @test A == B
    end
    for (key, _) in r_axis_writes(8)
        A = fresh_rbpot()
        ax = cc_ax(key...)
        apply_boundary_conditions_on_r_axis!(A, 1, ax, ax.interval, GBF)
        B = apply_op_pair!(fresh_rbpot(), 3, _ghost_ops(ax, 8, GBF))
        @test A == B
    end
end

@testset "r0 axis kernel matches the CPU broadcast implementation" begin
    # dim order (z, φ, r): 6 compressed z rows ⇒ up to 8 real z planes,
    # φ columns 1..5 (plus ghosts), the axis ring lives at r-row 2
    nnz = 8
    gw = rand(T, 6, 4)
    for even in (true, false), rbi in (1, 2)
        A = rand(T, 6, 7, 8, 2)
        B = copy(A)
        apply_boundary_conditions_at_r0!(A, rbi, nnz, gw, Val(even))
        backend = KernelAbstractions.CPU()
        r0_bc_kernel!(backend)(B, rbi, nnz, gw, even, ndrange = rbidx(nnz) - 1)
        KernelAbstractions.synchronize(backend)
        @test A ≈ B
    end
end

@testset "Geometrical axis weights array interface" begin
    w = rand(T, 4, 5)
    for W in (SolidStateDetectors.GeometricalCartesianAxisWeights{T},
              SolidStateDetectors.GeometricalAzimutalAxisWeights{T},
              SolidStateDetectors.GeometricalRadialAxisWeights{T})
        gw = W(w)
        @test size(gw) == size(w)
        @test gw[2, 3] == w[2, 3]
        @test all(gw[i] == w[i] for i in eachindex(w))
    end
end
