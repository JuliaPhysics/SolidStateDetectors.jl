

function apply_boundary_conditions_on_φ_axis!(  rbpot::Array{T, 4}, iz::Int, ir::Int, rbi::Int, nrbi::Int, ax::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :open, T},
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{only_2d})::Nothing where {T, only_2d}
    rbpot[iz,   1, ir, rbi] = rbpot[ iz, end - 1, ir, rbi] # cycling boundary
    rbpot[iz, end, ir, rbi] = rbpot[ iz,       2, ir, rbi] # cycling boundary
    nothing
end

function apply_boundary_conditions_on_φ_axis!(  rbpot::Array{T, 4}, iz::Int, ir::Int, rbi::Int, nrbi::Int, ax::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{false})::Nothing where {T}
    rbpot[iz,   1, ir, rbi] = rbpot[ iz,       3, ir, rbi] # cycling boundary
    rbpot[iz, end, ir, rbi] = rbpot[ iz, end - 2, ir, rbi] # cycling boundary
    nothing
end
function apply_boundary_conditions_on_φ_axis!(  rbpot::Array{T, 4}, iz::Int, ir::Int, rbi::Int, nrbi::Int, ax::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T}, only2d::Val{true})::Nothing where {T}
    rbpot[iz,   1, ir, nrbi] = rbpot[ iz, 2, ir, rbi] # cycling boundary
    rbpot[iz, end, ir, nrbi] = rbpot[ iz, 2, ir, rbi] # cycling boundary
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


function apply_boundary_conditions!(pssrb::PotentialSimulationSetupRB{T, 3, 4, Cylindrical, AT}, update_even_points::Val{even_points}, only2d::Val{only_2d}) where {T, even_points, only_2d, AT}
    rbi::Int = even_points ? rb_even::Int : rb_odd::Int
    nrbi::Int = even_points ? rb_odd::Int : rb_even::Int
    ax1 = pssrb.grid.axes[1]
    ax2 = pssrb.grid.axes[2]
    ax3 = pssrb.grid.axes[3]
    int1 = ax1.interval
    int2 = ax2.interval
    int3 = ax3.interval
    if only_2d
        iφ::Int = 2
        @inbounds for iz in axes(pssrb.potential, 1)
            # for ir in axes(pssrb.potential, 3)
            #     apply_boundary_conditions_on_φ_axis!( pssrb.potential, iz, ir, rbi, nrbi, ax2, int2, pssrb.grid_boundary_factors[2], only2d)
            # end
            apply_boundary_conditions_on_r_axis!( pssrb.potential, iz, iφ, rbi, ax1, int1, pssrb.grid_boundary_factors[1])
        end
        @inbounds for ir in axes(pssrb.potential, 3)
            apply_boundary_conditions_on_cyl_z_axis!( pssrb.potential, ir, iφ, rbi, ax3, int3, pssrb.grid_boundary_factors[3])
        end
    else
        @inbounds for iz in axes(pssrb.potential, 1)
            for iφ in axes(pssrb.potential, 2)
                apply_boundary_conditions_on_r_axis!( pssrb.potential, iz, iφ, rbi, ax1, int1, pssrb.grid_boundary_factors[1])
            end
            for ir in axes(pssrb.potential, 3)
                apply_boundary_conditions_on_φ_axis!( pssrb.potential, iz, ir, rbi, nrbi, ax2, int2, pssrb.grid_boundary_factors[2], only2d)
            end
        end
        @inbounds for iφ in axes(pssrb.potential, 2) # z boundaries
            for ir in axes(pssrb.potential, 3)
                apply_boundary_conditions_on_cyl_z_axis!( pssrb.potential, ir, iφ, rbi, ax3, int3, pssrb.grid_boundary_factors[3])
            end
        end
        begin # r = 0 handling
            nφ::Int = size(pssrb.potential, 2) - 1
            gw_φ::Array{T, 2} = pssrb.geom_weights[2].weights
            @inbounds for inz in 1:(size(pssrb.ϵ_r, 3) - 1)
                m::T = 0
                l::T = 0
                for inφ in 1:nφ
                    if even_points ? isodd(inz + inφ) : iseven(inz + inφ)
                        l += gw_φ[3, inφ] 
                        m += pssrb.potential[rbidx(inz), inφ + 1, 2, rbi] * gw_φ[3, inφ] 
                    end
                end
                m *= inv(l)
                for inφ in 1:nφ
                    if even_points ? isodd(inz + inφ) : iseven(inz + inφ)
                        pssrb.potential[rbidx(inz), inφ + 1, 2, rbi]::T = m
                    end
                end
            end
        end
    end
    nothing
end

