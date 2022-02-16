

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


function apply_boundary_conditions!(pssrb::PotentialCalculationSetup{T, Cylindrical}, update_even_points::Val{even_points}, only2d::Val{only_2d}) where {T, even_points, only_2d}
    rbi::Int = even_points ? rb_even::Int : rb_odd::Int
    # Radial Axis
    apply_boundary_conditions_on_r_axis!( pssrb.potential, rbi, pssrb.grid.axes[1], pssrb.grid.axes[1].interval, pssrb.grid_boundary_factors[1])
    # # Azimutal-Axis
    if !only_2d
        apply_boundary_conditions_on_φ_axis!( pssrb.potential, rbi, pssrb.grid.axes[2], pssrb.grid.axes[2].interval)
    end
    # Cylindrical Z-Axis -> same as Cartesian X-Axis
    apply_boundary_conditions_on_x_axis!( pssrb.potential, rbi, pssrb.grid.axes[3], pssrb.grid.axes[3].interval, pssrb.grid_boundary_factors[3])
    nothing
end

