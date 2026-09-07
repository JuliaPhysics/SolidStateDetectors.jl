function apply_boundary_conditions_on_z_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :infinite, :fixed}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[:, :, 1, rbi] .= grid_boundary_factors[1] .* view(rbpot, :, :, 2, rbi) # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_z_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :infinite, :infinite}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[:, :, 1, rbi] .= grid_boundary_factors[1] .* view(rbpot, :, :, 2, rbi)
    rbpot[:, :, end, rbi] .= grid_boundary_factors[2] .* view(rbpot, :, :, size(rbpot, 3) - 1, rbi)
    nothing
end
function apply_boundary_conditions_on_z_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :reflecting, :fixed}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[:, :, 1, rbi] .= view(rbpot, :, :, 3, rbi) # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_z_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :fixed, :infinite}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[:, :, end, rbi] .= grid_boundary_factors[2] .* view(rbpot, :, :, size(rbpot, 3) - 1, rbi) # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_z_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :fixed, :reflecting}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[:, :, end, rbi] .= view(rbpot, :, :, size(rbpot, 3) - 2, rbi) # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_z_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[:, :, 1, rbi] .= view(rbpot, :, :, 3, rbi) 
    rbpot[:, :, end, rbi] .= view(rbpot, :, :, size(rbpot, 3) - 2, rbi) 
    nothing
end
function apply_boundary_conditions_on_z_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :fixed, :fixed}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    nothing
end
function apply_boundary_conditions_on_z_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :open, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[:, :,   1, rbi] .= view(rbpot,  :, :, size(rbpot, 3) - 1, rbi) 
    rbpot[:, :, end, rbi] .= view(rbpot,  :, :,       2, rbi) 
    nothing
end


function apply_boundary_conditions_on_x_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :fixed, :fixed}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :infinite, :infinite}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, :, :, rbi] .= grid_boundary_factors[1] .* view(rbpot, 2, :, :, rbi)  # hmm, this is probably not fully correct since this is the red-black dimension
    rbpot[end, :, :, rbi] .= grid_boundary_factors[2] .* view(rbpot, size(rbpot, 1) - 1, :, :, rbi)  # anyhow its just an approximation
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :fixed, :infinite}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[end, :, :, rbi] .= grid_boundary_factors[2] .* view(rbpot, size(rbpot, 1) - 1, :, :, rbi)    # infinity boundaries in x (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :infinite, :fixed}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, :, :, rbi] .= grid_boundary_factors[1] .* view(rbpot, 2, :, :, rbi)  # infinity boundaries in x (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :infinite, :reflecting}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, :, :, rbi] .= grid_boundary_factors[1] .* view(rbpot, 2, :, :, rbi)  # infinity boundaries in x (only ±1 because this is the compressed (redblack) dimension)
    rbpot[end, :, :, rbi] .= view(rbpot, size(rbpot, 1) - 1, :, :, rbi) # (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :reflecting, :infinite}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[1, :, :, rbi]   .= view(rbpot, 2, :, :, rbi) #  (only ±1 because this is the compressed (redblack) dimension)
    rbpot[end, :, :, rbi] .= grid_boundary_factors[2] .* view(rbpot, size(rbpot, 1) - 1, :, :, rbi)    # infinity boundaries in x (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[1, :, :, rbi]   .= view(rbpot, 2, :, :, rbi) #  (only ±1 because this is the compressed (redblack) dimension)
    rbpot[end, :, :, rbi] .= view(rbpot, size(rbpot, 1) - 1, :, :, rbi) # (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :fixed, :reflecting}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[end, :, :, rbi] .= view(rbpot, size(rbpot, 1) - 1, :, :, rbi) # (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :reflecting, :fixed}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[1, :, :, rbi]   .= view(rbpot, 2, :, :, rbi) #  (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :open, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, :, :, rbi] .= view(rbpot, size(rbpot, 1), :, :, rbi) 
    rbpot[end, :, :, rbi] .= view(rbpot, 1, :, :, rbi) 
    nothing
end



function apply_boundary_conditions_on_y_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :fixed, :fixed}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :infinite, :infinite}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ :, 1, :, rbi]   .= grid_boundary_factors[1] .* view(rbpot, :, 2, :, rbi) 
    rbpot[ :, end, :, rbi] .= grid_boundary_factors[2] .* view(rbpot, :, size(rbpot, 2) - 1, :, rbi)   
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :fixed, :infinite}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ :, end, :, rbi] .= grid_boundary_factors[2] .* view(rbpot, :, size(rbpot, 2) - 1, :, rbi)   
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :infinite, :fixed}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ :, 1, :, rbi]   .= grid_boundary_factors[1] .* view(rbpot, :, 2, :, rbi) 
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ :, 1, :, rbi]   .= view(rbpot, :, 3, :, rbi) 
    rbpot[ :, end, :, rbi] .= view(rbpot, :, size(rbpot, 2) - 2, :, rbi)   
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :fixed, :reflecting}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ :, end, :, rbi] .= view(rbpot, :, size(rbpot, 2) - 2, :, rbi)   
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :reflecting, :fixed}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ :, 1, :, rbi]   .= view(rbpot, :, 3, :, rbi) 
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :infinite, :reflecting}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ :, 1, :, rbi]   .= grid_boundary_factors[1] .* view(rbpot, :, 2, :, rbi) 
    rbpot[ :, end, :, rbi] .= view(rbpot, :, size(rbpot, 2) - 2, :, rbi)   
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :reflecting, :infinite}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ :, 1, :, rbi]   .= view(rbpot, :, 3, :, rbi) 
    rbpot[ :, end, :, rbi] .= grid_boundary_factors[2] .* view(rbpot, :, size(rbpot, 2) - 1, :, rbi)   
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::AbstractArray{T, 4}, rbi::Int, ax::DiscreteAxis{T, :periodic, :periodic}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ :, 1, :, rbi]   .= view(rbpot, :, size(rbpot, 2) - 1, :, rbi) 
    rbpot[ :, end, :, rbi] .= view(rbpot, :, 2, :, rbi)   
    nothing
end


# Helper for apply_boundary_conditions! : shifts just the outer edge
# cells along `dim` (`margin` cells on each side, one red/black color) by
# ΔV, so the `:infinite` formula can be run in a gauge-shifted frame
# without touching the whole grid, this runs every SOR iteration, so a
# full-array shift would be a real slowdown.
#
# margin=3 matches the deepest read/write of any boundary variant in this
# file (`:reflecting` reads 3 cells in); the two edge ranges are kept from
# overlapping so a short axis isn't shifted twice.
@inline function _shift_axis_margin!(rbpot::AbstractArray{T, 4}, dim::Int, rbi::Int, ΔV::T; margin::Int = 3) where {T}
    n = size(rbpot, dim)
    lo_hi = min(margin, n)
    hi_lo = max(lo_hi + 1, n - margin + 1)
    if dim == 1
        @views rbpot[1:lo_hi, :, :, rbi] .+= ΔV
        hi_lo <= n && (@views rbpot[hi_lo:n, :, :, rbi] .+= ΔV)
    elseif dim == 2
        @views rbpot[:, 1:lo_hi, :, rbi] .+= ΔV
        hi_lo <= n && (@views rbpot[:, hi_lo:n, :, rbi] .+= ΔV)
    elseif dim == 3
        @views rbpot[:, :, 1:lo_hi, rbi] .+= ΔV
        hi_lo <= n && (@views rbpot[:, :, hi_lo:n, rbi] .+= ΔV)
    else
        error("_shift_axis_margin!: unsupported dim=$dim")
    end
    nothing
end

function apply_boundary_conditions!(pcs::PotentialCalculationSetup{T, Cartesian}, update_even_points::Val{even_points}, only2d::Val{only_2d}) where {T, even_points, only_2d}
    rbi::Int = even_points ? rb_even::Int : rb_odd::Int
    # `:infinite` decays each edge cell toward 0 (`factor * neighbor`),
    # this is correct only when V_ref is actually 0.
    # A gauge shift ({point=0,mantle=U} -> {point=-U,mantle=0}) moves V_ref to -U even though mantle
    # now reads 0. Fix: decay toward V_ref, not 0: n -> n - V_ref -> f*(n - V_ref) -> f*(n - V_ref) + V_ref
    V_ref = pcs.minimum_applied_potential
    _shift_axis_margin!(pcs.potential, 1, rbi, -V_ref)
    apply_boundary_conditions_on_x_axis!( pcs.potential, rbi, pcs.grid.axes[1], pcs.grid.axes[1].interval, pcs.grid_boundary_factors[1])
    _shift_axis_margin!(pcs.potential, 1, rbi, V_ref)
    _shift_axis_margin!(pcs.potential, 2, rbi, -V_ref)
    apply_boundary_conditions_on_y_axis!( pcs.potential, rbi, pcs.grid.axes[2], pcs.grid.axes[2].interval, pcs.grid_boundary_factors[2])
    _shift_axis_margin!(pcs.potential, 2, rbi, V_ref)
    _shift_axis_margin!(pcs.potential, 3, rbi, -V_ref)
    apply_boundary_conditions_on_z_axis!( pcs.potential, rbi, pcs.grid.axes[3], pcs.grid.axes[3].interval, pcs.grid_boundary_factors[3])
    _shift_axis_margin!(pcs.potential, 3, rbi, V_ref)
    nothing
end

