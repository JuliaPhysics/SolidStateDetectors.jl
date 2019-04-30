

function apply_boundary_conditions_on_φ_axis!(  rbpot::Array{T, 4}, iz::Int, ir::Int, rbi::Int, ax::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :open, T},
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{only_2d})::Nothing where {T, only_2d}
    rbpot[iz,   1, ir, rbi] = rbpot[ iz, end - 1, ir, rbi] # cycling boundary
    rbpot[iz, end, ir, rbi] = rbpot[ iz,       2, ir, rbi] # cycling boundary
    nothing
end

function apply_boundary_conditions_on_φ_axis!(  rbpot::Array{T, 4}, iz::Int, ir::Int, rbi::Int, ax::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{only_2d})::Nothing where {T, only_2d}
    if only_2d
        rbpot[iz,   1, ir, rbi] = rbpot[ iz, 2, ir, rbi] # cycling boundary
        rbpot[iz, end, ir, rbi] = rbpot[ iz, 2, ir, rbi] # cycling boundary
    else
        rbpot[iz,   1, ir, rbi] = rbpot[ iz, 3, ir, rbi] # cycling boundary
        rbpot[iz, end, ir, rbi] = rbpot[ iz, end - 2, ir, rbi] # cycling boundary
    end
    nothing
end
function apply_boundary_conditions_on_r_axis!(  rbpot::Array{T, 4}, iz::Int, iφ::Int, rbi::Int, ax::DiscreteAxis{T, :r0, :infinite}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[iz, iφ, end, rbi] = grid_boundary_factors[2] * rbpot[iz, iφ, end - 2, rbi] # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_r_axis!(  rbpot::Array{T, 4}, iz::Int, iφ::Int, rbi::Int, ax::DiscreteAxis{T, :r0, :reflecting}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[iz, iφ, end, rbi] = rbpot[iz, iφ, end - 2, rbi] # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_r_axis!(  rbpot::Array{T, 4}, iz::Int, iφ::Int, rbi::Int, ax::DiscreteAxis{T, :r0, :fixed}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    # rbpot[iz, iφ, end, rbi] = rbpot[iz, iφ, end - 2, rbi] # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_cyl_z_axis!(  rbpot::Array{T, 4}, ir::Int, iφ::Int, rbi::Int, ax::DiscreteAxis{T, :infinite, :infinite}, int::Interval{:closed, :closed, T},
                                                    grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, iφ, ir, rbi] = grid_boundary_factors[1] * rbpot[ 2, iφ, ir, rbi]  # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    rbpot[end, iφ, ir, rbi] = grid_boundary_factors[2] * rbpot[ end - 1, iφ, ir, rbi]    # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_cyl_z_axis!(  rbpot::Array{T, 4}, ir::Int, iφ::Int, rbi::Int, ax::DiscreteAxis{T, :infinite, :fixed}, int::Interval{:closed, :closed, T},
                                                    grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, iφ, ir, rbi] = grid_boundary_factors[1] * rbpot[ 2, iφ, ir, rbi]  # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_cyl_z_axis!(  rbpot::Array{T, 4}, ir::Int, iφ::Int, rbi::Int, ax::DiscreteAxis{T, :fixed, :infinite}, int::Interval{:closed, :closed, T},
                                                    grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[end, iφ, ir, rbi] = grid_boundary_factors[2] * rbpot[ end - 1, iφ, ir, rbi]    # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_cyl_z_axis!(  rbpot::Array{T, 4}, ir::Int, iφ::Int, rbi::Int, ax::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval{:closed, :closed, T},
                                                    grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, iφ, ir, rbi] = rbpot[ 2, iφ, ir, rbi]  # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    rbpot[end, iφ, ir, rbi] = rbpot[ end - 1, iφ, ir, rbi]    # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_cyl_z_axis!(  rbpot::Array{T, 4}, ir::Int, iφ::Int, rbi::Int, ax::DiscreteAxis{T, :fixed, :reflecting}, int::Interval{:closed, :closed, T},
                                                    grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[end, iφ, ir, rbi] = rbpot[ end - 1, iφ, ir, rbi]    # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_cyl_z_axis!(  rbpot::Array{T, 4}, ir::Int, iφ::Int, rbi::Int, ax::DiscreteAxis{T, :reflecting, :fixed}, int::Interval{:closed, :closed, T},
                                                    grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, iφ, ir, rbi] = rbpot[ 2, iφ, ir, rbi]  # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_cyl_z_axis!(  rbpot::Array{T, 4}, ir::Int, iφ::Int, rbi::Int, ax::DiscreteAxis{T, :fixed, :fixed}, int::Interval{:closed, :closed, T},
                                                    grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    nothing
end


function apply_boundary_conditions!(fssrb::PotentialSimulationSetupRB{T, 3, 4, :cylindrical}, update_even_points::Val{even_points}, only2d::Val{only_2d}) where {T, even_points, only_2d}
    rbi::Int = even_points ? rb_even::Int : rb_odd::Int
    if only_2d
        iφ::Int = 2
        @inbounds for iz in axes(fssrb.potential, 1)
            for ir in axes(fssrb.potential, 3)
                apply_boundary_conditions_on_φ_axis!( fssrb.potential, iz, ir, rbi, fssrb.grid.axes[2], fssrb.grid.axes[2].interval, fssrb.grid_boundary_factors[2], only2d)
            end
            apply_boundary_conditions_on_r_axis!( fssrb.potential, iz, iφ, rbi, fssrb.grid.axes[1], fssrb.grid.axes[1].interval, fssrb.grid_boundary_factors[1])
        end
        @inbounds for ir in axes(fssrb.potential, 3)
            apply_boundary_conditions_on_cyl_z_axis!( fssrb.potential, ir, iφ, rbi, fssrb.grid.axes[3], fssrb.grid.axes[3].interval, fssrb.grid_boundary_factors[3])
        end
    else
        @inbounds for iz in axes(fssrb.potential, 1)
            for iφ in axes(fssrb.potential, 2)
                apply_boundary_conditions_on_r_axis!( fssrb.potential, iz, iφ, rbi, fssrb.grid.axes[1], fssrb.grid.axes[1].interval, fssrb.grid_boundary_factors[1])
            end
            for ir in axes(fssrb.potential, 3)
                apply_boundary_conditions_on_φ_axis!( fssrb.potential, iz, ir, rbi, fssrb.grid.axes[2], fssrb.grid.axes[2].interval, fssrb.grid_boundary_factors[2], only2d)
            end
        end
        @inbounds for iφ in axes(fssrb.potential, 2) # z boundaries
            for ir in axes(fssrb.potential, 3)
                apply_boundary_conditions_on_cyl_z_axis!( fssrb.potential, ir, iφ, rbi, fssrb.grid.axes[3], fssrb.grid.axes[3].interval, fssrb.grid_boundary_factors[3])
            end
        end
        begin # r = 0 handling
            nφ::Int = size(fssrb.potential, 2) - 1
            gw_φ::Array{T, 2} = fssrb.geom_weights[2].weights
            @inbounds for inz in 1:(size(fssrb.ϵ, 3) - 1)
                m::T = 0
                l::T = 0
                for inφ in 1:nφ
                    if even_points ? isodd(inz + inφ) : iseven(inz + inφ)
                        l += gw_φ[3, inφ] 
                        m += fssrb.potential[rbidx(inz), inφ + 1, 2, rbi] * gw_φ[3, inφ] 
                    end
                end
                m *= inv(l)
                for inφ in 1:nφ
                    if even_points ? isodd(inz + inφ) : iseven(inz + inφ)
                        fssrb.potential[rbidx(inz), inφ + 1, 2, rbi]::T = m
                    end
                end
            end
        end
    end
    nothing
end

