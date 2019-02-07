function apply_boundary_conditions_on_θ_axis!(  rbpot::Array{T, 4}, iz::Int, ir::Int, rbi::Int, ax::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :open, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[iz,   1, ir, rbi] = rbpot[ iz, end - 1, ir, rbi] # cycling boundary
    rbpot[iz, end, ir, rbi] = rbpot[ iz,       2, ir, rbi] # cycling boundary
    nothing
end
function apply_boundary_conditions_on_θ_axis!(  rbpot::Array{T, 4}, iz::Int, ir::Int, rbi::Int, ax::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[iz,   1, ir, rbi] = rbpot[ iz, end - 1, ir, rbi] # cycling boundary
    rbpot[iz, end, ir, rbi] = rbpot[ iz,       2, ir, rbi] # cycling boundary
    nothing
end
function apply_boundary_conditions_on_r_axis!(  rbpot::Array{T, 4}, iz::Int, iθ::Int, rbi::Int, ax::DiscreteAxis{T, :r0, :infinite}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[iz, iθ, end, rbi] = grid_boundary_factors[2] * rbpot[iz, iθ, end - 2, rbi] # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_r_axis!(  rbpot::Array{T, 4}, iz::Int, iθ::Int, rbi::Int, ax::DiscreteAxis{T, :r0, :reflecting}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[iz, iθ, end, rbi] = rbpot[iz, iθ, end - 2, rbi] # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_cyl_z_axis!(  rbpot::Array{T, 4}, ir::Int, iθ::Int, rbi::Int, ax::DiscreteAxis{T, :infinite, :infinite}, int::Interval{:closed, :closed, T},
                                                    grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, iθ, ir, rbi] = grid_boundary_factors[1] * rbpot[ 2, iθ, ir, rbi]  # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    rbpot[end, iθ, ir, rbi] = grid_boundary_factors[2] * rbpot[ end - 1, iθ, ir, rbi]    # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_cyl_z_axis!(  rbpot::Array{T, 4}, ir::Int, iθ::Int, rbi::Int, ax::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval{:closed, :closed, T},
                                                    grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, iθ, ir, rbi] = rbpot[ 2, iθ, ir, rbi]  # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    rbpot[end, iθ, ir, rbi] = rbpot[ end - 1, iθ, ir, rbi]    # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    nothing
end


function apply_boundary_conditions!(fssrb::PotentialSimulationSetupRB{T, 3, 4, :Cylindrical}, update_even_points::Val{even_points}, only2d::Val{only_2d}) where {T, even_points, only_2d}
    rbi::Int = even_points ? rb_even::Int : rb_odd::Int
    if only_2d
        iθ::Int = 2
        @inbounds for iz in axes(fssrb.potential, 1)
            for ir in axes(fssrb.potential, 3)
                apply_boundary_conditions_on_θ_axis!( fssrb.potential, iz, ir, rbi, fssrb.grid.axes[2], fssrb.grid.axes[2].interval, fssrb.grid_boundary_factors[2])
            end
            apply_boundary_conditions_on_r_axis!( fssrb.potential, iz, iθ, rbi, fssrb.grid.axes[1], fssrb.grid.axes[1].interval, fssrb.grid_boundary_factors[1])
        end
        @inbounds for ir in axes(fssrb.potential, 3)
            apply_boundary_conditions_on_cyl_z_axis!( fssrb.potential, ir, iθ, rbi, fssrb.grid.axes[3], fssrb.grid.axes[3].interval, fssrb.grid_boundary_factors[3])
        end
    else
        @inbounds for iz in axes(fssrb.potential, 1)
            for iθ in axes(fssrb.potential, 2)
                apply_boundary_conditions_on_r_axis!( fssrb.potential, iz, iθ, rbi, fssrb.grid.axes[1], fssrb.grid.axes[1].interval, fssrb.grid_boundary_factors[1])
            end
            for ir in axes(fssrb.potential, 3)
                apply_boundary_conditions_on_θ_axis!( fssrb.potential, iz, ir, rbi, fssrb.grid.axes[2], fssrb.grid.axes[2].interval, fssrb.grid_boundary_factors[2])
            end
        end
        @inbounds for iθ in axes(fssrb.potential, 2) # z boundaries
            for ir in axes(fssrb.potential, 3)
                apply_boundary_conditions_on_cyl_z_axis!( fssrb.potential, ir, iθ, rbi, fssrb.grid.axes[3], fssrb.grid.axes[3].interval, fssrb.grid_boundary_factors[3])
            end
        end
        begin # r = 0 handling
            nθ::Int = size(fssrb.potential, 2) - 1
            gw_θ::Array{T, 2} = fssrb.geom_weights[2].weights
            @inbounds for inz in 1:(size(fssrb.ϵ, 3) - 1)
                m::T = 0
                l::T = 0
                for inθ in 1:nθ
                    if even_points ? isodd(inz + inθ) : iseven(inz + inθ)
                        l += gw_θ[3, inθ] 
                        m += fssrb.potential[rbidx(inz), inθ + 1, 2, rbi] * gw_θ[3, inθ] 
                    end
                end
                m *= inv(l)
                for inθ in 1:nθ
                    if even_points ? isodd(inz + inθ) : iseven(inz + inθ)
                        fssrb.potential[rbidx(inz), inθ + 1, 2, rbi]::T = m
                    end
                end
            end
        end
    end
    nothing
end

