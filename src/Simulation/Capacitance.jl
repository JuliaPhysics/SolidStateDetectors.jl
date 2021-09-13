export calculate_mutual_capacitance

function calculate_mutual_capacitance(sim::Simulation, ij::Tuple{Int, Int})
    calculate_mutual_capacitance(sim.weighting_potentials[ij[1]], sim.weighting_potentials[ij[2]], sim.ϵ_r)
end

function calculate_mutual_capacitance(
            wpi::WeightingPotential{T,3,CS}, wpj::WeightingPotential{T,3,CS},
            ϵ_r::DielectricDistribution{T,3,CS};
            consider_multiplicity::Bool = true ) where {T, CS}
    c_ij::T = zero(T)

    int_p1 = interpolated_scalarfield(wpi)
    int_p2 = interpolated_scalarfield(wpj)
    int_ϵ_r = interpolated_scalarfield(ϵ_r)
    cylindrical = CS == Cylindrical
    phi_2D = cylindrical && size(wpi, 2) == size(wpj, 2) == 1
    grid = phi_2D ? get_2π_potential(wpi, n_points_in_φ = 2).grid : _get_closed_potential(wpi).grid
    grid_mps = get_extended_midpoints_grid(grid)
    for i3 in 1:size(grid, 3)-1
        for i2 in 1:size(grid, 2)-1
            for i1 in 1:size(grid, 1)-1
                w1, w2, w3 = voxel_widths(grid, i1, i2, i3)
                dV = voxel_volume(grid, i1, i2, i3, w1, w2, w3)

                pt_voxel_mid = GridPoint(grid_mps, (i1 + 1, i2 + 1, i3 + 1))
                ϵ_r_voxel = get_interpolation(int_ϵ_r, pt_voxel_mid, CS)

                efs = (Vector{T}(undef, 3), Vector{T}(undef, 3))
                for (i, int_p) in enumerate((int_p1, int_p2))
                    p000 = get_interpolation(int_p, GridPoint(grid, (i1    , i2    , i3    )), CS)
                    p100 = get_interpolation(int_p, GridPoint(grid, (i1 + 1, i2    , i3    )), CS)
                    p010 = get_interpolation(int_p, GridPoint(grid, (i1    , i2 + 1, i3    )), CS)
                    p110 = get_interpolation(int_p, GridPoint(grid, (i1 + 1, i2 + 1, i3    )), CS)
                    p001 = get_interpolation(int_p, GridPoint(grid, (i1    , i2    , i3 + 1)), CS)
                    p101 = get_interpolation(int_p, GridPoint(grid, (i1 + 1, i2    , i3 + 1)), CS)
                    p011 = get_interpolation(int_p, GridPoint(grid, (i1    , i2 + 1, i3 + 1)), CS)
                    p111 = get_interpolation(int_p, GridPoint(grid, (i1 + 1, i2 + 1, i3 + 1)), CS)
                    efv1 = ( (p100 - p000) + (p110 - p010) + (p101 - p001) + (p111 - p011) ) / (4 * w1)
                    efv2 = if cylindrical
                        _w2 = (grid[2].ticks[i2 + 1] - grid[2].ticks[i2])
                        if i1 == 1
                            ((p110 - p100)/(_w2*grid[1].ticks[i1+1]) +
                             (p111 - p101)/(_w2*grid[1].ticks[i1+1]) ) / 2
                        else
                            ((p010 - p000)/(_w2*grid[1].ticks[i1]) +
                             (p110 - p100)/(_w2*grid[1].ticks[i1+1]) +
                             (p011 - p001)/(_w2*grid[1].ticks[i1]) +
                             (p111 - p101)/(_w2*grid[1].ticks[i1+1])) / 4
                        end
                    else
                        ( (p010 - p000) + (p110 - p100) + (p011 - p001) + (p111 - p101) ) / (4 * w2)
                    end
                    efv3 = ( (p001 - p000) + (p101 - p100) + (p011 - p010) + (p111 - p110) ) / (4 * w3)
                    efs[i][:] = [efv1, efv2, efv3]
                end
                c_ij += sum(efs[1] .* efs[2]) * dV * ϵ_r_voxel
            end
        end
    end
    phi_2D && (c_ij *= 2)
    consider_multiplicity && (c_ij *= multiplicity(grid))
    return uconvert(u"pF", c_ij*u"m" * ϵ0*u"F/m" )
end

export calculate_capacitance_matrix
"""
    calculate_capacitance_matrix(sim::Simulation{T}) where {T}

Calculates the Maxwell Capacitance Matrix.
"""
function calculate_capacitance_matrix(sim::Simulation{T}) where {T}
    @assert !ismissing(sim.weighting_potentials) "The weighting_potentials needs to be calculated first."
    n = length(sim.weighting_potentials)
    C = zeros(typeof(one(T) * u"pF"), (n, n))
    for i in 1:n
        for j in i:n
            C[j, i] = C[i, j] = calculate_mutual_capacitance(sim, (i, j))
        end
    end
    return C
end