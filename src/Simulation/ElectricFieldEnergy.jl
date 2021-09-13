export calculate_stored_energy
"""
    calculate_stored_energy(sim::Simulation; consider_multiplicity::Bool = true)

Calculates and returns the energy stored in the [`ElectricField`](@ref) of a 
[`SolidStateDetector`](@ref) in a given [`Simulation`](@ref) in units of J.

## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) with `sim.detector` for which the stored energy is calculated.

## Keywords 
* `consider_multiplicity::Bool = true`: Whether symmetries of the system should be taken into account. 
    For example, in case of true coaxial detector center around the origin and calculated on a cartesian grid 
    with the `x-axis` going from `[0, x_max]` and the `y-axis` going from `[0, y_max]` the multiplicity is 4
    and, if `consider_multiplicity == true`, the returned value is already multiplied by 4.
"""
function calculate_stored_energy(sim::Simulation; consider_multiplicity::Bool = true)
    @assert !ismissing(sim.electric_potential) "Electric potential has not been calculated yet. Please run `calculate_electric_potential!(sim)` first."
    calculate_stored_energy(sim.electric_potential, sim.ϵ_r; consider_multiplicity)
end


function calculate_stored_energy(ep::ElectricPotential{T,3,CS}, ϵ::DielectricDistribution{T,3,CS}; consider_multiplicity::Bool = true) where {T <: SSDFloat, CS}
    cylindrical = CS == Cylindrical
    phi_2D = cylindrical && size(ep, 2) == 1
    ep3d = phi_2D ? get_2π_potential(ep, n_points_in_φ = 2) : _get_closed_potential(ep)
    grid = ep3d.grid
    W::T = 0
    for i3 in 1:size(grid, 3)-1
        for i2 in 1:size(grid, 2)-1
            for i1 in 1:size(grid, 1)-1
                w1, w2, w3 = voxel_widths(grid, i1, i2, i3)
                dV = voxel_volume(grid, i1, i2, i3, w1, w2, w3)

                _ϵ = ϵ.data[i1 + 1, i2 + 1, i3 + 1]

                ep000 = ep3d.data[i1    , i2    , i3    ]
                ep100 = ep3d.data[i1 + 1, i2    , i3    ]
                ep010 = ep3d.data[i1    , i2 + 1, i3    ]
                ep110 = ep3d.data[i1 + 1, i2 + 1, i3    ]
                ep001 = ep3d.data[i1    , i2    , i3 + 1]
                ep101 = ep3d.data[i1 + 1, i2    , i3 + 1]
                ep011 = ep3d.data[i1    , i2 + 1, i3 + 1]
                ep111 = ep3d.data[i1 + 1, i2 + 1, i3 + 1]

                efv1 = ( (ep100 - ep000) + (ep110 - ep010) + (ep101 - ep001) + (ep111 - ep011) ) / (4 * w1)
                efv2 = if cylindrical
                    _w2 = (grid[2].ticks[i2 + 1] - grid[2].ticks[i2])
                    if i1 == 1
                        ((ep110 - ep100)/(_w2*grid[1].ticks[i1+1]) +
                        (ep111 - ep101)/(_w2*grid[1].ticks[i1+1])) / 2
                    else
                        ((ep010 - ep000)/(_w2*grid[1].ticks[i1]) +
                         (ep110 - ep100)/(_w2*grid[1].ticks[i1+1]) +
                         (ep011 - ep001)/(_w2*grid[1].ticks[i1]) +
                         (ep111 - ep101)/(_w2*grid[1].ticks[i1+1])) / 4
                    end
                else
                    ( (ep010 - ep000) + (ep110 - ep100) + (ep011 - ep001) + (ep111 - ep101) ) / (4 * w2)
                end
                efv3 = ( (ep001 - ep000) + (ep101 - ep100) + (ep011 - ep010) + (ep111 - ep110) ) / (4 * w3)
                W += sum((efv1, efv2, efv3).^2) * dV * _ϵ
            end
        end
    end
    E = W * ϵ0 / 2 * u"J"
    phi_2D && (E *= 2)
    return consider_multiplicity ? E * multiplicity(grid) : E
end
