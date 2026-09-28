

function apply_boundary_conditions_on_φ_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, 
                                                ax::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :open, T} )::Nothing where {T}
    @inbounds rbpot[:,   1, :, rbi] .= view(rbpot, :, size(rbpot, 2) - 1, :, rbi) # cycling boundary
    @inbounds rbpot[:, end, :, rbi] .= view(rbpot, :,                  2, :, rbi) # cycling boundary
    nothing
end

function apply_boundary_conditions_on_φ_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, 
                                                ax::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval{:closed, :closed, T} )::Nothing where {T}
    @inbounds rbpot[:,   1, :, rbi] .= view(rbpot, :,                  3, :, rbi) # cycling boundary
    @inbounds rbpot[:, end, :, rbi] .= view(rbpot, :, size(rbpot, 2) - 2, :, rbi) # cycling boundary
    nothing
end


function apply_boundary_conditions_on_r_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :r0, :infinite}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    @inbounds rbpot[:, :, end, rbi] .= grid_boundary_factors[2] .* view(rbpot, :, :, size(rbpot, 3) - 2, rbi) # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_r_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :r0, :reflecting}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    @inbounds rbpot[:, :, end, rbi] .= view(rbpot, :, :, size(rbpot, 3) - 2, rbi) # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_r_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :r0, :fixed}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    nothing
end

# One work-item per z-row: sum the row's φ columns of matching parity with the
# precomputed weights and broadcast the result back to those columns. Replaces
# a getindex + gemv + broadcast sequence (a temporary and several kernel
# launches per SOR iteration) on GPU backends.
@kernel function r0_bc_kernel!(rbpot, rbi::Int, nnz::Int, gw, even_points::Bool)
    irbz = @index(Global, Linear)
    rbz = irbz + 1
    inφ_ub = size(rbpot, 2) - 1
    @inbounds begin
        last_odd = nnz - (iseven(nnz) ? 1 : 0)   # odd real-z rows
        if rbz <= div(last_odd, 2) + mod(last_odd, 2) + 1
            gwrow = even_points ? 6 : 5
            φ0 = even_points ? 2 : 1
            s = zero(eltype(rbpot))
            j = 0
            for inφ in φ0:2:inφ_ub
                j += 1
                s += rbpot[rbz, inφ + 1, 2, rbi] * gw[gwrow, j]
            end
            for inφ in φ0:2:inφ_ub
                rbpot[rbz, inφ + 1, 2, rbi] = s
            end
        end
        last_even = nnz - (isodd(nnz) ? 1 : 0)   # even real-z rows
        if rbz <= div(last_even, 2) + mod(last_even, 2) + 1
            gwrow = even_points ? 5 : 6
            φ0 = even_points ? 1 : 2
            s = zero(eltype(rbpot))
            j = 0
            for inφ in φ0:2:inφ_ub
                j += 1
                s += rbpot[rbz, inφ + 1, 2, rbi] * gw[gwrow, j]
            end
            for inφ in φ0:2:inφ_ub
                rbpot[rbz, inφ + 1, 2, rbi] = s
            end
        end
    end
end

function apply_boundary_conditions_at_r0!(
    rbpot::GPUArrays.AbstractGPUArray{T, 4},
    rbi::Int,
    nnz::Int,
    gw::AbstractArray{T},
    update_even_points::Val{even_points}
    ) where {T, even_points}
    kern = r0_bc_kernel!(_ka_get_backend(rbpot))
    kern(rbpot, rbi, nnz, gw, even_points, ndrange = rbidx(nnz) - 1)
    nothing
end

function apply_boundary_conditions_at_r0!(
    rbpot::AbstractArray{T, 4},
    rbi::Int,
    nnz::Int,
    gw::AbstractArray{T},
    update_even_points::Val{even_points}
    ) where {T, even_points}
    if even_points
        inz_range = 1:2:nnz # inz odd
        rbz_range = 2:rbidx(inz_range.stop)
        inφ_range = 2:2:size(rbpot, 2)-1 # must be even 
        inφ_range_pot = inφ_range .+ 1
        @inbounds @views rbpot[rbz_range, inφ_range_pot, 2, rbi] .= rbpot[rbz_range, inφ_range_pot, 2, rbi] * gw[6, 1:length(inφ_range)]
        inz_range = 2:2:nnz # inz even
        rbz_range = 2:rbidx(inz_range.stop)
        inφ_range = 1:2:size(rbpot, 2)-1 # must be odd 
        inφ_range_pot = inφ_range .+ 1
        @inbounds @views rbpot[rbz_range, inφ_range_pot, 2, rbi] .= rbpot[rbz_range, inφ_range_pot, 2, rbi] * gw[5, 1:length(inφ_range)]
    else
        inz_range = 1:2:nnz # inz odd
        rbz_range = 2:rbidx(inz_range.stop)
        inφ_range = 1:2:size(rbpot, 2)-1 # must be odd 
        inφ_range_pot = inφ_range .+ 1
        @inbounds @views rbpot[rbz_range, inφ_range_pot, 2, rbi] .= rbpot[rbz_range, inφ_range_pot, 2, rbi] * gw[5, 1:length(inφ_range)]
        inz_range = 2:2:nnz # inz even
        rbz_range = 2:rbidx(inz_range.stop)
        inφ_range = 2:2:size(rbpot, 2)-1 # must be even 
        inφ_range_pot = inφ_range .+ 1
        @inbounds @views rbpot[rbz_range, inφ_range_pot, 2, rbi] .= rbpot[rbz_range, inφ_range_pot, 2, rbi] * gw[6, 1:length(inφ_range)]
    end
    nothing
end

function apply_boundary_conditions!(pcs::PotentialCalculationSetup{T, Cylindrical}, update_even_points::Val{even_points}, only2d::Val{only_2d}) where {T, even_points, only_2d}
    rbi::Int = even_points ? rb_even::Int : rb_odd::Int
    # Radial Axis at rMax
    apply_boundary_conditions_on_r_axis!( pcs.potential, rbi, pcs.grid.axes[1], pcs.grid.axes[1].interval, pcs.grid_boundary_factors[1])
    if !only_2d
        # Radial axis at r0
        apply_boundary_conditions_at_r0!(pcs.potential, rbi, size(pcs.ϵ_r, 3)-1, pcs.geom_weights[2], update_even_points)
        # Azimutal-Axis
        apply_boundary_conditions_on_φ_axis!( pcs.potential, rbi, pcs.grid.axes[2], pcs.grid.axes[2].interval)
    end
    # Cylindrical Z-Axis -> same as Cartesian X-Axis
    apply_boundary_conditions_on_x_axis!( pcs.potential, rbi, pcs.grid.axes[3], pcs.grid.axes[3].interval, pcs.grid_boundary_factors[3])
    nothing
end

