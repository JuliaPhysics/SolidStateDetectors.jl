abstract type AbstractPotentialSimulationSetup{T, N} end

# """
#     PotentialSimulationSetup{T, N, S} <: AbstractPotentialSimulationSetup{T, N}
# 
# Collection struct. It holds the grid, the potential, the point types, the charge density and the dielectric distribution.
# """
struct PotentialSimulationSetup{T, N, S} <: AbstractPotentialSimulationSetup{T, N}
    grid::Grid{T, N, S}
    potential::Array{T, N}
    pointtypes::Array{PointType, N}
    q_eff_imp::Array{T, N}
    q_eff_fix::Array{T, N}
    ϵ_r::Array{T, N}
end

struct PotentialSimulationSetupRB{T, N1, N2, S, TGW, AT} <: AbstractPotentialSimulationSetup{T, N1}
    grid::Grid{T, N1, S, AT}
    potential::Array{T, N2}
    pointtypes::Array{PointType, N2}
    volume_weights::Array{T, N2}
    q_eff_imp::Array{T, N2}
    q_eff_fix::Array{T, N2}
    ϵ_r::Array{T, N1}
    geom_weights::TGW   
    sor_const::Array{T, 1}
    bias_voltage::T
    maximum_applied_potential::T
    minimum_applied_potential::T
    depletion_handling_potential_limit::T
    grid_boundary_factors::NTuple{3, NTuple{2, T}}
end

function sizeof(fssrb::PotentialSimulationSetup{T, N})::Int where {T, N}
    s::Int = sizeof(fssrb.grid)
    s += sizeof(fssrb.pointtypes)
    s += sizeof(fssrb.potential)
    s += sizeof(fssrb.ϵ_r)
    s += sizeof(fssrb.q_eff_imp)
    s += sizeof(fssrb.q_eff_fix)
    return s
end

function sizeof(fssrb::PotentialSimulationSetupRB{T, N1, N2})::Int where {T, N1, N2}
    s::Int = sizeof(fssrb.grid)
    s += sizeof(fssrb.pointtypes)
    s += sizeof(fssrb.potential)
    s += sizeof(fssrb.volume_weights)
    s += sizeof(fssrb.ϵ_r)
    s += sizeof(fssrb.q_eff_imp)
    s += sizeof(fssrb.q_eff_fix)
    for idim in 1:N1
        s += sizeof(fssrb.geom_weights[idim].weights)
    end
    s += sizeof(fssrb.sor_const)
    s += sizeof(fssrb.bias_voltage)
    s += sizeof(fssrb.maximum_applied_potential)
    s += sizeof(fssrb.minimum_applied_potential)
    s += sizeof(fssrb.depletion_handling_potential_limit)
    s += sizeof(fssrb.grid_boundary_factors)
    return s
end

include("BoundaryConditions/BoundaryConditions.jl")
include("PotentialSimulationSetupRBCylindrical.jl")
include("PotentialSimulationSetupRBCartesian3D.jl")
