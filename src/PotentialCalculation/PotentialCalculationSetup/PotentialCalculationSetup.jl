abstract type AbstractPotentialCalculationSetup{T, N} end

"""
    struct PotentialCalculationSetup

`PotentialCalculationSetup` holds the grid, fields and certain precalculated fixed parameters for the field calculation.
This struct will be calculated after each refinement as is depends on the grid. 

`grid`: The 3-dimensional grid, either Cartesian or cylindrical, on which the field will be calculated. 

The fields:
* `potential`: This is the 4-dimensional array of the extended 3D potential array. 
The fourth dimensions comes from the red-black (even-odd) division in order to parallelize the field calculation. 
Extended means that the grid holds one additional tick at both sides of all axes necessary for boundary handling (reflecting, fixed, ...) at the ends of the grid.
* `point_types`: Also a 4-dimensional array. Same structure as the `potential`-array.
* `volume_weights`: Also a 4-dimensional array. Same structure as the `potential`-array.
* `q_eff_imp`: Also a 4-dimensional array. Same structure as the `potential`-array.
* `q_eff_fix`: Also a 4-dimensional array. Same structure as the `potential`-array.
* `ϵ_r`: Normal 3-dimensional array! In order to calculate the final 6 weights for each grid point in the SOR, 
the eight values of the dielectric permittivity (of each octant) around the grid point needs to be loaded. 
Here not special division is possible for the red-black (even-odd) division. 

Precalculated parameters:
* `geom_weights`: The parts of the calculation of the six weights in the SOR can be precalculated. Those are stored here for each axis/dimension.
* `sor_const`: Vector holding the SOR constants. In the cartesian case only the first entry is used. As the optimal value for the SOR constant
depends on the grid, the constant is linear increased and the array holds the respective value for each radial axis tick. 
* `bias_voltage`: `maximum_applied_potential - minimum_applied_potential`. Used for depletion handling, but might be obsolete by now. 
* `maximum_applied_potential`: Used for depletion handling, but might be obsolete by now. 
* `minimum_applied_potential`: Used for depletion handling, but might be obsolete by now. 
* `grid_boundary_factors`: Used in the application of boundary conditions in the field calculation for decaying (infinite) boundary conditions 
to approximate the decay of the potential (depending on the grid).
"""
struct PotentialCalculationSetup{
            T, S, N1, DATN1<:AbstractArray{T,N1}, 
            N2, AT, 
            DATN2<:AbstractArray{T,N2}, 
            DATPT<:AbstractArray{PointType,N2}, 
            DATGW<:AbstractArray{T,2},
            DATSOR<:AbstractArray{T,1}
        } <: AbstractPotentialCalculationSetup{T, N1}
    grid::Grid{T, N1, S, AT}
    potential::DATN2 # Array{T, N2} or e.g. CuArray{T, N2}
    point_types::DATPT # Array{PointType, N2} or ...
    volume_weights::DATN2 # Array{T, N2}
    q_eff_imp::DATN2 # Array{T, N2}
    q_eff_fix::DATN2 # Array{T, N2}
    ϵ_r::DATN1 # Array{T, N1}
    geom_weights::NTuple{3, DATGW} # NTuple{3, Array{T, 2}}
    sor_const::DATSOR # Array{T, 1}
    bias_voltage::T
    maximum_applied_potential::T
    minimum_applied_potential::T    
    grid_boundary_factors::NTuple{3, NTuple{2, T}}
end

function Adapt.adapt_structure(to, pcs::PotentialCalculationSetup{T, S, 3}) where {T, S}
    PotentialCalculationSetup(
        pcs.grid,
        adapt(to, pcs.potential),
        adapt(to, pcs.point_types),
        adapt(to, pcs.volume_weights),
        adapt(to, pcs.q_eff_imp),
        adapt(to, pcs.q_eff_fix),
        adapt(to, pcs.ϵ_r),
        adapt(to, pcs.geom_weights),
        adapt(to, pcs.sor_const),
        pcs.bias_voltage,
        pcs.maximum_applied_potential,
        pcs.minimum_applied_potential,
        pcs.grid_boundary_factors
    )
end

include("BoundaryConditions/BoundaryConditions.jl")
include("PotentialCalculationSetupCylindrical.jl")
include("PotentialCalculationSetupCartesian3D.jl")
