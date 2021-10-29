@doc raw"""
    calculate_mutual_capacitance(sim::Simulation, ij::Tuple{Int, Int}; consider_multiplicity::Bool = true)

Returns the mutual capacitance between the contacts with ID `i = ij[1]` and `j = ij[2]`.
It is calculated via the weighting potentials of the contacts, ``\Phi_i^w(\vec{r})`` and ``\Phi_j^w(\vec{r})``:
```math
c_{ij} = \epsilon_0 \int_{World} \nabla \Phi_i^w(\vec{r}) ϵ_r(\vec{r}) \nabla \Phi_j^w(\vec{r}) d\vec{r}
```

!!! note
    These are elements of the Mawell Capcitance Matrix. Look up [Capacitances](@ref) for more information.

!!! note 
    The electric potential as well as the two weighting potentials of both contacts have to be calculated.

## Arguments
* `sim::Simulation`: [`Simulation`](@ref) for which the capacitance matrix is calculated.
* `ij::Tuple{Int,Int}`: Tuple of indices of the contacts for which the capacitance should be calculated. 

## Keywords 
* `consider_multiplicity::Bool = true`: Whether symmetries of the system should be taken into account. 
    For example, in case of true coaxial detector center around the origin and calculated on a cartesian grid 
    with the `x-axis` going from `[0, x_max]` and the `y-axis` going from `[0, y_max]` the multiplicity is 4
    and, if `consider_multiplicity == true`, the returned value is already multiplied by 4.
"""
function calculate_mutual_capacitance(sim::Simulation, ij::Tuple{Int, Int}; consider_multiplicity::Bool = true)
    _calculate_mutual_capacitance(sim.weighting_potentials[ij[1]], sim.weighting_potentials[ij[2]], sim.ϵ_r; consider_multiplicity)
end


function _calculate_mutual_capacitance(
            wpi::WeightingPotential{T,3,CS}, wpj::WeightingPotential{T,3,CS},
            ϵ_r::DielectricDistribution{T,3,CS};
            consider_multiplicity::Bool = true ) where {T, CS}

    int_p1 = interpolated_scalarfield(wpi)
    int_p2 = interpolated_scalarfield(wpj)
    int_ϵ_r = interpolated_scalarfield(ϵ_r)
    cylindrical = CS == Cylindrical
    phi_2D = cylindrical && size(wpi, 2) == size(wpj, 2) == 1
    grid = phi_2D ? get_2π_potential(wpi, n_points_in_φ = 2).grid : _get_closed_potential(wpi).grid
    grid_mps = get_extended_midpoints_grid(grid)
    
    c_ij::T = _calculate_mutual_capacitance(grid, grid_mps, int_ϵ_r, int_p1, int_p2)

    phi_2D && (c_ij *= 2)
    consider_multiplicity && (c_ij *= multiplicity(grid))
    return uconvert(u"pF", c_ij*u"m" * ϵ0*u"F/m" )
end

function _calculate_mutual_capacitance(grid::Grid{T, 3, CS}, grid_mps, int_ϵ_r, int_p1, int_p2) where {T, CS}
    c_ij::T = zero(T)
    for i3 in 1:size(grid, 3)-1
        for i2 in 1:size(grid, 2)-1
            for i1 in 1:size(grid, 1)-1
                w1, w2, w3 = voxel_widths(grid, i1, i2, i3)
                dV = voxel_volume(grid, i1, i2, i3, w1, w2, w3)

                pt_voxel_mid = GridPoint(grid_mps, (i1 + 1, i2 + 1, i3 + 1))
                ϵ_r_voxel = get_interpolation(int_ϵ_r, pt_voxel_mid, CS)

                efs_1 = _approximate_potential_gradient(int_p1, grid, i1, i2, i3, w1, w2, w3) 
                efs_2 = _approximate_potential_gradient(int_p2, grid, i1, i2, i3, w1, w2, w3) 

                c_ij += sum(efs_1 .* efs_2) * dV * ϵ_r_voxel
            end
        end
    end
    return c_ij
end

function _approximate_potential_gradient(int_p, grid::Grid{T, 3, CS}, i1, i2, i3, w1, w2, w3) where {T, CS}
    p000 = get_interpolation(int_p, GridPoint(grid, (i1    , i2    , i3    )), CS)
    p100 = get_interpolation(int_p, GridPoint(grid, (i1 + 1, i2    , i3    )), CS)
    p010 = get_interpolation(int_p, GridPoint(grid, (i1    , i2 + 1, i3    )), CS)
    p110 = get_interpolation(int_p, GridPoint(grid, (i1 + 1, i2 + 1, i3    )), CS)
    p001 = get_interpolation(int_p, GridPoint(grid, (i1    , i2    , i3 + 1)), CS)
    p101 = get_interpolation(int_p, GridPoint(grid, (i1 + 1, i2    , i3 + 1)), CS)
    p011 = get_interpolation(int_p, GridPoint(grid, (i1    , i2 + 1, i3 + 1)), CS)
    p111 = get_interpolation(int_p, GridPoint(grid, (i1 + 1, i2 + 1, i3 + 1)), CS)
    efv1 = ( (p100 - p000) + (p110 - p010) + (p101 - p001) + (p111 - p011) ) / (4 * w1)
    efv2 = if CS == Cylindrical
        _w2 = (grid.axes[2].ticks[i2 + 1] - grid.axes[2].ticks[i2])
        if i1 == 1
            ((p110 - p100)/(_w2*grid.axes[1].ticks[i1+1]) +
                (p111 - p101)/(_w2*grid.axes[1].ticks[i1+1]) ) / 2
        else
            ((p010 - p000)/(_w2*grid.axes[1].ticks[i1]) +
                (p110 - p100)/(_w2*grid.axes[1].ticks[i1+1]) +
                (p011 - p001)/(_w2*grid.axes[1].ticks[i1]) +
                (p111 - p101)/(_w2*grid.axes[1].ticks[i1+1])) / 4
        end
    else
        ( (p010 - p000) + (p110 - p100) + (p011 - p001) + (p111 - p101) ) / (4 * w2)
    end
    efv3 = ( (p001 - p000) + (p101 - p100) + (p011 - p010) + (p111 - p110) ) / (4 * w3)
    return (efv1, efv2, efv3)
end



@doc raw"""
    calculate_capacitance_matrix(sim::Simulation{T}; consider_multiplicity::Bool = true) where {T}

Calculates the Maxwell Capacitance `N×N`-Matrix in units of pF,
where `N` is the number of contacts of `sim.detector`.
The individual elements, ``c_{i,j}``, are calculated via 
[`calculate_mutual_capacitance(sim::Simulation, (i,j)::Tuple{Int,Int})`](@ref).
The matrix should be symmetric. The difference of `C[i,j]` and `C[j,i]` are due 
to numerical precision in the integration due to the different grids of the two weighting potentials.

## Arguments
* `sim::Simulation`: [`Simulation`](@ref) for which the capacitance matrix is calculated.

## Keywords 
* `consider_multiplicity::Bool = true`: Whether symmetries of the system should be taken into account. 
    For example, in case of true coaxial detector center around the origin and calculated on a cartesian grid 
    with the `x-axis` going from `[0, x_max]` and the `y-axis` going from `[0, y_max]` the multiplicity is 4
    and, if `consider_multiplicity == true`, the returned value is already multiplied by 4.
"""
function calculate_capacitance_matrix(sim::Simulation{T}; consider_multiplicity::Bool = true) where {T}
    @assert !ismissing(sim.ϵ_r) "The electric potential needs to be calculated first."
    @assert !ismissing(sim.weighting_potentials) "The weighting_potentials needs to be calculated first."
    n = length(sim.weighting_potentials)
    C = zeros(typeof(one(T) * u"pF"), (n, n))
    for i in 1:n
        for j in 1:n
            C[j, i] = if !ismissing(sim.weighting_potentials[i]) && !ismissing(sim.weighting_potentials[j]) 
                calculate_mutual_capacitance(sim, (i, j); consider_multiplicity)
            else
                missing
            end
        end
    end
    return C
end