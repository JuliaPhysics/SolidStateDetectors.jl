function apply_boundary_conditions_on_z_axis!(  rbpot::Array{T, 4}, ix::Int, iy::Int, rbi::Int, ax::DiscreteAxis{T, :infinite, :fixed}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ix, iy, 1, rbi] = grid_boundary_factors[1] * rbpot[ix, iy, 2, rbi] # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_z_axis!(  rbpot::Array{T, 4}, ix::Int, iy::Int, rbi::Int, ax::DiscreteAxis{T, :infinite, :infinite}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ix, iy, 1, rbi] = grid_boundary_factors[1] * rbpot[ix, iy, 2, rbi] 
    rbpot[ix, iy, end, rbi] = grid_boundary_factors[2] * rbpot[ix, iy, end - 1, rbi]
    nothing
end
function apply_boundary_conditions_on_z_axis!(  rbpot::Array{T, 4}, ix::Int, iy::Int, rbi::Int, ax::DiscreteAxis{T, :reflecting, :fixed}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ix, iy, 1, rbi] = rbpot[ix, iy, 3, rbi] # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_z_axis!(  rbpot::Array{T, 4}, ix::Int, iy::Int, rbi::Int, ax::DiscreteAxis{T, :fixed, :infinite}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ix, iy, end, rbi] = grid_boundary_factors[2] * rbpot[ix, iy, end - 1, rbi] # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_z_axis!(  rbpot::Array{T, 4}, ix::Int, iy::Int, rbi::Int, ax::DiscreteAxis{T, :fixed, :reflecting}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ix, iy, end, rbi] = rbpot[ix, iy, end - 2, rbi] # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_z_axis!(  rbpot::Array{T, 4}, ix::Int, iy::Int, rbi::Int, ax::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ix, iy, 1, rbi] = rbpot[ix, iy, 3, rbi] 
    rbpot[ix, iy, end, rbi] = rbpot[ix, iy, end - 2, rbi] 
    nothing
end
function apply_boundary_conditions_on_z_axis!(  rbpot::Array{T, 4}, ix::Int, iy::Int, rbi::Int, ax::DiscreteAxis{T, :fixed, :fixed}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    nothing
end
function apply_boundary_conditions_on_z_axis!(  rbpot::Array{T, 4}, ix::Int, iy::Int, rbi::Int, ax::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :open, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[ix, iy,   1, rbi] = rbpot[ ix, iy, end - 1, rbi] 
    rbpot[ix, iy, end, rbi] = rbpot[ ix, iy,       2, rbi] 
    nothing
end


function apply_boundary_conditions_on_x_axis!(  rbpot::Array{T, 4}, iy::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :fixed, :fixed}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::Array{T, 4}, iy::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :infinite, :infinite}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, iy, iz, rbi] = grid_boundary_factors[1] * rbpot[ 2, iy, iz, rbi]  # hmm, this is probably not fully correct since this is the red-black dimension
    rbpot[end, iy, iz, rbi] = grid_boundary_factors[2] * rbpot[ end - 1, iy, iz, rbi]  # anyhow its just an approximation
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::Array{T, 4}, iy::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :fixed, :infinite}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[end, iy, iz, rbi] = grid_boundary_factors[2] * rbpot[ end - 1, iy, iz, rbi]    # infinity boundaries in x (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::Array{T, 4}, iy::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :infinite, :fixed}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, iy, iz, rbi] = grid_boundary_factors[1] * rbpot[ 2, iy, iz, rbi]  # infinity boundaries in x (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::Array{T, 4}, iy::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :infinite, :reflecting}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, iy, iz, rbi] = grid_boundary_factors[1] * rbpot[ 2, iy, iz, rbi]  # infinity boundaries in x (only ±1 because this is the compressed (redblack) dimension)
    rbpot[end, iy, iz, rbi] = rbpot[end - 1, iy, iz, rbi] # (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::Array{T, 4}, iy::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :reflecting, :infinite}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[1, iy, iz, rbi]   = rbpot[2, iy, iz, rbi] #  (only ±1 because this is the compressed (redblack) dimension)
    rbpot[end, iy, iz, rbi] = grid_boundary_factors[2] * rbpot[ end - 1, iy, iz, rbi]    # infinity boundaries in x (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::Array{T, 4}, iy::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[1, iy, iz, rbi]   = rbpot[2, iy, iz, rbi] #  (only ±1 because this is the compressed (redblack) dimension)
    rbpot[end, iy, iz, rbi] = rbpot[end - 1, iy, iz, rbi] # (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::Array{T, 4}, iy::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :fixed, :reflecting}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[end, iy, iz, rbi] = rbpot[end - 1, iy, iz, rbi] # (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::Array{T, 4}, iy::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :reflecting, :fixed}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[1, iy, iz, rbi]   = rbpot[2, iy, iz, rbi] #  (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_x_axis!(  rbpot::Array{T, 4}, iy::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :open, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, iy, iz, rbi] = rbpot[end, iy, iz, rbi] 
    rbpot[end, iy, iz, rbi] = rbpot[  1, iy, iz, rbi] 
    nothing
end



function apply_boundary_conditions_on_y_axis!(  rbpot::Array{T, 4}, ix::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :fixed, :fixed}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{only_2d})::Nothing where {T, only_2d}
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::Array{T, 4}, ix::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :infinite, :infinite}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{only_2d})::Nothing where {T, only_2d}
    rbpot[ ix, 1, iz, rbi]   = grid_boundary_factors[1] * rbpot[ ix, 2, iz, rbi] 
    rbpot[ ix, end, iz, rbi] = grid_boundary_factors[2] * rbpot[ ix, end - 1, iz, rbi]   
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::Array{T, 4}, ix::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :fixed, :infinite}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{only_2d})::Nothing where {T, only_2d}
    rbpot[ ix, end, iz, rbi] = grid_boundary_factors[2] * rbpot[ ix, end - 1, iz, rbi]   
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::Array{T, 4}, ix::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :infinite, :fixed}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{only_2d})::Nothing where {T, only_2d}
    rbpot[ ix, 1, iz, rbi]   = grid_boundary_factors[1] * rbpot[ ix, 2, iz, rbi] 
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::Array{T, 4}, ix::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{only_2d})::Nothing where {T, only_2d}
    rbpot[ ix, 1, iz, rbi]   = rbpot[ ix, 3, iz, rbi] 
    rbpot[ ix, end, iz, rbi] = rbpot[ ix, end - 2, iz, rbi]   
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::Array{T, 4}, ix::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :fixed, :reflecting}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{only_2d})::Nothing where {T, only_2d}
    rbpot[ ix, end, iz, rbi] = rbpot[ ix, end - 2, iz, rbi]   
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::Array{T, 4}, ix::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :reflecting, :fixed}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{only_2d})::Nothing where {T, only_2d}
    rbpot[ ix, 1, iz, rbi]   = rbpot[ ix, 3, iz, rbi] 
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::Array{T, 4}, ix::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :infinite, :reflecting}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{only_2d})::Nothing where {T, only_2d}
    rbpot[ ix, 1, iz, rbi]   = grid_boundary_factors[1] * rbpot[ ix, 2, iz, rbi] 
    rbpot[ ix, end, iz, rbi] = rbpot[ ix, end - 2, iz, rbi]   
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::Array{T, 4}, ix::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :reflecting, :infinite}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{only_2d})::Nothing where {T, only_2d}
    rbpot[ ix, 1, iz, rbi]   = rbpot[ ix, 3, iz, rbi] 
    rbpot[ ix, end, iz, rbi] = grid_boundary_factors[2] * rbpot[ ix, end - 1, iz, rbi]   
    nothing
end
function apply_boundary_conditions_on_y_axis!(  rbpot::Array{T, 4}, ix::Int, iz::Int, rbi::Int, ax::DiscreteAxis{T, :periodic, :periodic}, int::Interval,
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{only_2d})::Nothing where {T, only_2d}
    rbpot[ ix, 1, iz, rbi]   = rbpot[ ix, end - 1, iz, rbi] 
    rbpot[ ix, end, iz, rbi] = rbpot[ ix, 2, iz, rbi]   
    nothing
end

function apply_boundary_conditions!(fssrb::PotentialSimulationSetupRB{T, N1, N2, :cartesian}, update_even_points::Val{even_points}, only2d::Val{only_2d}) where {T, N1, N2, even_points, only_2d}
    rbi::Int = even_points ? rb_even::Int : rb_odd::Int

    if only_2d
        iy::Int = 2
        error("Boundary handling for 2D simulation (x&z) in cartesian coordinated not yet implemented.")
    else
        @inbounds for ix in axes(fssrb.potential, 1)
            for iy in axes(fssrb.potential, 2)
                apply_boundary_conditions_on_z_axis!( fssrb.potential, ix, iy, rbi, fssrb.grid.axes[3], fssrb.grid.axes[3].interval, fssrb.grid_boundary_factors[3])
            end
            for iz in axes(fssrb.potential, 3)
                apply_boundary_conditions_on_y_axis!( fssrb.potential, ix, iz, rbi, fssrb.grid.axes[2], fssrb.grid.axes[2].interval, fssrb.grid_boundary_factors[2], only2d)
            end
        end
        @inbounds for iy in axes(fssrb.potential, 2) # z boundaries
            for iz in axes(fssrb.potential, 3)
                apply_boundary_conditions_on_x_axis!( fssrb.potential, iy, iz, rbi, fssrb.grid.axes[1], fssrb.grid.axes[1].interval, fssrb.grid_boundary_factors[1])
            end
        end
    end
    nothing
end


