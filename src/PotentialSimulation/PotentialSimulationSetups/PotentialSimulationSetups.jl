abstract type AbstractPotentialSimulationSetup{T, N} end

struct PotentialSimulationSetup{T, N, S} <: AbstractPotentialSimulationSetup{T, N}
    grid::Grid{T, N, S}
    potential::Array{T, N}
    pointtypes::Array{PointType, N}
    ρ::Array{T, N}
    ϵ::Array{T, N}
end

struct PotentialSimulationSetupRB{T, N1, N2, S} <: AbstractPotentialSimulationSetup{T, N1}
    grid::Grid{T, N1, S}
    potential::Array{T, N2}
    pointtypes::Array{PointType, N2}
    volume_weights::Array{T, N2}
    ρ::Array{T, N2}
    ϵ::Array{T, N1}
    geom_weights::NTuple{N1, AbstractGeometricalAxisWeights{T}}        
    sor_const::Array{T, 1}
    bias_voltage::T
    maximum_applied_potential::T
    minimum_applied_potential::T
    depletion_handling_potential_limit::T
    bulk_is_ptype::Bool
    grid_boundary_factors::NTuple{3, NTuple{2, T}}
end

function sizeof(fssrb::PotentialSimulationSetup{T, N})::Int where {T, N}
    s::Int = sizeof(fssrb.grid)
    s += sizeof(fssrb.pointtypes)
    s += sizeof(fssrb.potential)
    s += sizeof(fssrb.ϵ)
    s += sizeof(fssrb.ρ)
    return s
end

function sizeof(fssrb::PotentialSimulationSetupRB{T, N1, N2})::Int where {T, N1, N2}
    s::Int = sizeof(fssrb.grid)
    s += sizeof(fssrb.pointtypes)
    s += sizeof(fssrb.potential)
    s += sizeof(fssrb.volume_weights)
    s += sizeof(fssrb.ϵ)
    s += sizeof(fssrb.ρ)
    for idim in 1:N1
        s += sizeof(fssrb.geom_weights[idim].weights)
    end
    s += sizeof(fssrb.sor_const)
    s += sizeof(fssrb.bias_voltage)
    s += sizeof(fssrb.maximum_applied_potential)
    s += sizeof(fssrb.minimum_applied_potential)
    s += sizeof(fssrb.depletion_handling_potential_limit)
    s += sizeof(fssrb.bulk_is_ptype)
    s += sizeof(fssrb.grid_boundary_factors)
    return s
end

function PotentialSimulationSetup(ssd::SolidStateDetector{T}, grid::Grid{T, N, S} = Grid(ssd), potential_array::Union{Missing, Array{T, N}} = missing)::PotentialSimulationSetup{T, N, S} where {T, N, S}   
    fssrb::PotentialSimulationSetupRB{T, N, N + 1, S} = PotentialSimulationSetupRB(ssd, grid, potential_array)
    return PotentialSimulationSetup{T, N, S}( Grid(fssrb), ElectricPotentialArray(fssrb), PointTypeArray(fssrb), ChargeDensityArray(fssrb), DielektrikumDistributionArray(fssrb)  )
end

include("PotentialSimulationSetupRBCylindrical.jl")


@recipe function f( pss::PotentialSimulationSetup{T, 3, :Cylindrical};
                    r = missing,
                    θ = missing,
                    z = missing,
                    n_points_in_θ = 36 ) where {T}
    g::Grid{T, 3, :Cylindrical} = pss.grid
    layout --> (2, 2) 

    cross_section::Symbol, idx::Int = if ismissing(θ) && ismissing(r) && ismissing(z)
        :θ, 1
    elseif !ismissing(θ) && ismissing(r) && ismissing(z)
        θ_rad::T = T(deg2rad(θ))
        while !(g[:θ].interval.left <= θ_rad <= g[:θ].interval.right)
            if θ_rad > g[:θ].interval.right
                θ_rad -= g[:θ].interval.right - g[:θ].interval.left
            elseif θ_rad < g[:θ].interval.left
                θ_rad += g[:θ].interval.right - g[:θ].interval.left
            end
        end
        :θ, searchsortednearest(g[:θ], θ_rad)
    elseif ismissing(θ) && !ismissing(r) && ismissing(z)
        :r, searchsortednearest(g[:r], T(r))
    elseif ismissing(θ) && ismissing(r) && !ismissing(z)
        :z, searchsortednearest(g[:z], T(z))
    else
        error(ArgumentError, ": Only one of the keywords `r, θ, z` is allowed.")
    end
    value::T = if cross_section == :θ
        g[:θ][idx]
    elseif cross_section == :r    
        g[:r][idx]
    elseif cross_section == :z
        g[:z][idx]
    end

    if cross_section == :θ
        θ --> θ
    elseif cross_section == :z
        z --> z
    elseif cross_section == :r
        r --> r
    end

    @series begin
        subplot := 1
        PointTypes(pss)
    end
    @series begin
        subplot := 2
        ChargeDensity(pss)
    end
    @series begin
        subplot := 3
        DielectricDistribution(pss)
    end
    @series begin
        subplot := 4
        ElectricPotential(pss, n_points_in_θ=n_points_in_θ)
    end
end
