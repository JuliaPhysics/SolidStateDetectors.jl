# GPU boundary-condition path: the CPU implementation applies the ghost-face
# updates as ~10 view broadcasts per red-black color, whose per-call host
# overhead dominates the GPU iteration loop. Here each face update is encoded
# as (active, source index, factor) with ghost = factor * pot[source] and all
# six faces are applied by one kernel launch per color.
#
# The op tables below are 1:1 transcriptions of the broadcast methods in
# BoundaryConditionsCylindrical.jl / BoundaryConditionsCartesian.jl.
# Ghost entries with more than one ghost coordinate are never read by the SOR
# stencil, so the faces may be applied concurrently.

const _GhostFaceOp{T} = Tuple{Bool, Int, T}
_no_ghost_op(::Type{T}) where {T} = (false, 1, one(T))

# faces of the compressed (red-black) dimension: neighbors sit at rows 2 / end-1
_ghost_ops_compressed(::DiscreteAxis{T, :fixed, :fixed}, n::Int, gbf) where {T} =
    (_no_ghost_op(T), _no_ghost_op(T))
_ghost_ops_compressed(::DiscreteAxis{T, :infinite, :infinite}, n::Int, gbf) where {T} =
    ((true, 2, gbf[1]), (true, n - 1, gbf[2]))
_ghost_ops_compressed(::DiscreteAxis{T, :fixed, :infinite}, n::Int, gbf) where {T} =
    (_no_ghost_op(T), (true, n - 1, gbf[2]))
_ghost_ops_compressed(::DiscreteAxis{T, :infinite, :fixed}, n::Int, gbf) where {T} =
    ((true, 2, gbf[1]), _no_ghost_op(T))
_ghost_ops_compressed(::DiscreteAxis{T, :infinite, :reflecting}, n::Int, gbf) where {T} =
    ((true, 2, gbf[1]), (true, n - 1, one(T)))
_ghost_ops_compressed(::DiscreteAxis{T, :reflecting, :infinite}, n::Int, gbf) where {T} =
    ((true, 2, one(T)), (true, n - 1, gbf[2]))
_ghost_ops_compressed(::DiscreteAxis{T, :reflecting, :reflecting}, n::Int, gbf) where {T} =
    ((true, 2, one(T)), (true, n - 1, one(T)))
_ghost_ops_compressed(::DiscreteAxis{T, :fixed, :reflecting}, n::Int, gbf) where {T} =
    (_no_ghost_op(T), (true, n - 1, one(T)))
_ghost_ops_compressed(::DiscreteAxis{T, :reflecting, :fixed}, n::Int, gbf) where {T} =
    ((true, 2, one(T)), _no_ghost_op(T))
_ghost_ops_compressed(::DiscreteAxis{T, :periodic, :periodic}, n::Int, gbf) where {T} =
    ((true, n - 1, one(T)), (true, 2, one(T)))

# faces of the non-compressed dimensions: reflecting neighbors sit at rows 3 / end-2
_ghost_ops(::DiscreteAxis{T, :fixed, :fixed}, n::Int, gbf) where {T} =
    (_no_ghost_op(T), _no_ghost_op(T))
_ghost_ops(::DiscreteAxis{T, :infinite, :infinite}, n::Int, gbf) where {T} =
    ((true, 2, gbf[1]), (true, n - 1, gbf[2]))
_ghost_ops(::DiscreteAxis{T, :fixed, :infinite}, n::Int, gbf) where {T} =
    (_no_ghost_op(T), (true, n - 1, gbf[2]))
_ghost_ops(::DiscreteAxis{T, :infinite, :fixed}, n::Int, gbf) where {T} =
    ((true, 2, gbf[1]), _no_ghost_op(T))
_ghost_ops(::DiscreteAxis{T, :infinite, :reflecting}, n::Int, gbf) where {T} =
    ((true, 2, gbf[1]), (true, n - 2, one(T)))
_ghost_ops(::DiscreteAxis{T, :reflecting, :infinite}, n::Int, gbf) where {T} =
    ((true, 3, one(T)), (true, n - 1, gbf[2]))
_ghost_ops(::DiscreteAxis{T, :reflecting, :reflecting}, n::Int, gbf) where {T} =
    ((true, 3, one(T)), (true, n - 2, one(T)))
_ghost_ops(::DiscreteAxis{T, :fixed, :reflecting}, n::Int, gbf) where {T} =
    (_no_ghost_op(T), (true, n - 2, one(T)))
_ghost_ops(::DiscreteAxis{T, :reflecting, :fixed}, n::Int, gbf) where {T} =
    ((true, 3, one(T)), _no_ghost_op(T))
_ghost_ops(::DiscreteAxis{T, :periodic, :periodic}, n::Int, gbf) where {T} =
    ((true, n - 1, one(T)), (true, 2, one(T)))
# the radial axis: no left face (r = 0 is handled separately), and the
# infinite/reflecting right face reads row end-2 (as in the broadcast methods)
_ghost_ops(::DiscreteAxis{T, :r0, :infinite}, n::Int, gbf) where {T} =
    (_no_ghost_op(T), (true, n - 2, gbf[2]))
_ghost_ops(::DiscreteAxis{T, :r0, :reflecting}, n::Int, gbf) where {T} =
    (_no_ghost_op(T), (true, n - 2, one(T)))
_ghost_ops(::DiscreteAxis{T, :r0, :fixed}, n::Int, gbf) where {T} =
    (_no_ghost_op(T), _no_ghost_op(T))

@kernel function faces_bc_kernel!(rbpot, rbi::Int,
        a1l::Bool, s1l::Int, f1l, a1r::Bool, s1r::Int, f1r,
        a2l::Bool, s2l::Int, f2l, a2r::Bool, s2r::Int, f2r,
        a3l::Bool, s3l::Int, f3l, a3r::Bool, s3r::Int, f3r)
    i, j = @index(Global, NTuple)
    n1 = size(rbpot, 1); n2 = size(rbpot, 2); n3 = size(rbpot, 3)
    @inbounds begin
        if i <= n1 && j <= n2
            if a3l rbpot[i, j,  1, rbi] = f3l * rbpot[i, j, s3l, rbi] end
            if a3r rbpot[i, j, n3, rbi] = f3r * rbpot[i, j, s3r, rbi] end
        end
        if i <= n1 && j <= n3
            if a2l rbpot[i,  1, j, rbi] = f2l * rbpot[i, s2l, j, rbi] end
            if a2r rbpot[i, n2, j, rbi] = f2r * rbpot[i, s2r, j, rbi] end
        end
        if i <= n2 && j <= n3
            if a1l rbpot[ 1, i, j, rbi] = f1l * rbpot[s1l, i, j, rbi] end
            if a1r rbpot[n1, i, j, rbi] = f1r * rbpot[s1r, i, j, rbi] end
        end
    end
end

function _apply_faces_bc!(rbpot::GPUArrays.AbstractGPUArray{T, 4}, rbi::Int,
        ops1::NTuple{2, _GhostFaceOp{T}}, ops2::NTuple{2, _GhostFaceOp{T}}, ops3::NTuple{2, _GhostFaceOp{T}}) where {T}
    kern = faces_bc_kernel!(_ka_get_backend(rbpot))
    n1, n2, n3 = size(rbpot, 1), size(rbpot, 2), size(rbpot, 3)
    kern(rbpot, rbi,
        ops1[1]..., ops1[2]..., ops2[1]..., ops2[2]..., ops3[1]..., ops3[2]...,
        ndrange = (max(n1, n2), max(n2, n3)))
    nothing
end

function apply_boundary_conditions!(pcs::PotentialCalculationSetup{T, Cylindrical, 3, <:GPUArrays.AbstractGPUArray},
        update_even_points::Val{even_points}, only2d::Val{only_2d}) where {T, even_points, only_2d}
    rbi::Int = even_points ? rb_even::Int : rb_odd::Int
    rbpot = pcs.potential
    if !only_2d
        apply_boundary_conditions_at_r0!(rbpot, rbi, size(pcs.ϵ_r, 3) - 1, pcs.geom_weights[2], update_even_points)
    end
    # red-black dimension order is (z, φ, r)
    ops_z = _ghost_ops_compressed(pcs.grid.axes[3], size(rbpot, 1), pcs.grid_boundary_factors[3])
    ops_φ = only_2d ? (_no_ghost_op(T), _no_ghost_op(T)) : _ghost_ops(pcs.grid.axes[2], size(rbpot, 2), pcs.grid_boundary_factors[2])
    ops_r = _ghost_ops(pcs.grid.axes[1], size(rbpot, 3), pcs.grid_boundary_factors[1])
    _apply_faces_bc!(rbpot, rbi, ops_z, ops_φ, ops_r)
    nothing
end

function apply_boundary_conditions!(pcs::PotentialCalculationSetup{T, Cartesian, 3, <:GPUArrays.AbstractGPUArray},
        update_even_points::Val{even_points}, only2d::Val{only_2d}) where {T, even_points, only_2d}
    rbi::Int = even_points ? rb_even::Int : rb_odd::Int
    rbpot = pcs.potential
    ops_x = _ghost_ops_compressed(pcs.grid.axes[1], size(rbpot, 1), pcs.grid_boundary_factors[1])
    ops_y = _ghost_ops(pcs.grid.axes[2], size(rbpot, 2), pcs.grid_boundary_factors[2])
    ops_z = _ghost_ops(pcs.grid.axes[3], size(rbpot, 3), pcs.grid_boundary_factors[3])
    _apply_faces_bc!(rbpot, rbi, ops_x, ops_y, ops_z)
    nothing
end
