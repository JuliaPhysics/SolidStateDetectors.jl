abstract type AbstractPotentialSimulationSetup{T, N} end

struct PotentialSimulationSetupRB{
            T, S, N1, DATN1<:AbstractArray{T,N1}, 
            N2, AT, 
            DATN2<:AbstractArray{T,N2}, 
            DATPT<:AbstractArray{PointType,N2}, 
            DATGW<:AbstractArray{T,2},
            DATSOR<:AbstractArray{T,1}
        } <: AbstractPotentialSimulationSetup{T, N1}
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
    depletion_handling_potential_limit::T
    grid_boundary_factors::NTuple{3, NTuple{2, T}}
end

function Adapt.adapt_structure(to, pssrb::PotentialSimulationSetupRB{T, S, 3}) where {T, S}
    PotentialSimulationSetupRB(
        pssrb.grid,
        adapt(to, pssrb.potential),
        adapt(to, pssrb.point_types),
        adapt(to, pssrb.volume_weights),
        adapt(to, pssrb.q_eff_imp),
        adapt(to, pssrb.q_eff_fix),
        adapt(to, pssrb.ϵ_r),
        adapt(to, pssrb.geom_weights),
        adapt(to, pssrb.sor_const),
        pssrb.bias_voltage,
        pssrb.maximum_applied_potential,
        pssrb.minimum_applied_potential,
        pssrb.depletion_handling_potential_limit,
        pssrb.grid_boundary_factors
    )
end

include("BoundaryConditions/BoundaryConditions.jl")
include("PotentialSimulationSetupRBCylindrical.jl")
include("PotentialSimulationSetupRBCartesian3D.jl")
