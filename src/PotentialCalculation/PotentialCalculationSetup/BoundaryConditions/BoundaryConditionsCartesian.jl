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


function apply_boundary_conditions!(pcs::PotentialCalculationSetup{T, Cartesian}, update_even_points::Val{even_points}, only2d::Val{only_2d}) where {T, even_points, only_2d}
    rbi::Int = even_points ? rb_even::Int : rb_odd::Int
    apply_boundary_conditions_on_x_axis!( pcs.potential, rbi, pcs.grid.axes[1], pcs.grid.axes[1].interval, pcs.grid_boundary_factors[1])
    apply_boundary_conditions_on_y_axis!( pcs.potential, rbi, pcs.grid.axes[2], pcs.grid.axes[2].interval, pcs.grid_boundary_factors[2])
    apply_boundary_conditions_on_z_axis!( pcs.potential, rbi, pcs.grid.axes[3], pcs.grid.axes[3].interval, pcs.grid_boundary_factors[3])
    nothing
end

