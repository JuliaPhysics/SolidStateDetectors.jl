struct PotentialSimulationSetupRBGPU{T, N1, N2, S, TGW, AT} <: AbstractPotentialSimulationSetup{T, N1}
    grid::Grid{T, N1, S, AT}
    potential::CuArray{T, N2}
    point_types::CuArray{PointType, N2}
    volume_weights::CuArray{T, N2}
    q_eff_imp::CuArray{T, N2}
    q_eff_fix::CuArray{T, N2}
    ϵ_r::CuArray{T, N1}
    geom_weights::NTuple{3, CuArray{T, 2}}
    sor_const::CuArray{T, 1}
    bias_voltage::T
    maximum_applied_potential::T
    minimum_applied_potential::T
    depletion_handling_potential_limit::T
    grid_boundary_factors::NTuple{3, NTuple{2, T}}
end

function ElectricPotentialArray(pssrb::PotentialSimulationSetupRBGPU{T, 3, 4, Cylindrical})::Array{T, 3} where {T}
    pot::Array{T, 3} = Array{T, 3}(undef, size(pssrb.grid))
    rbpot = Array(pssrb.potential);
    for iz in axes(pot, 3)
        irbz::Int = rbidx(iz)
        for iφ in axes(pot, 2)
            irbφ::Int = iφ + 1
            idxsum::Int = iz + iφ
            for ir in axes(pot, 1)
                irbr::Int = ir + 1
                rbi::Int = iseven(idxsum + ir) ? rb_even::Int : rb_odd::Int
                pot[ir, iφ, iz] = rbpot[ irbz, irbφ, irbr, rbi ]
            end
        end
    end
    return pot
end

function PotentialSimulationSetupRBGPU(pssrb::PotentialSimulationSetupRB{T, N1, N2, S, TGW, AT}) where {T, N1, N2, S, TGW, AT}
    PotentialSimulationSetupRBGPU{T, N1, N2, S, TGW, AT}(
        pssrb.grid,
        CuArray(pssrb.potential),
        CuArray(pssrb.point_types),
        CuArray(pssrb.volume_weights),
        CuArray(pssrb.q_eff_imp),
        CuArray(pssrb.q_eff_fix),
        CuArray(pssrb.ϵ_r),
        broadcast(w -> CuArray(w), pssrb.geom_weights),
        CuArray(pssrb.sor_const),
        pssrb.bias_voltage,
        pssrb.maximum_applied_potential,
        pssrb.minimum_applied_potential,
        pssrb.depletion_handling_potential_limit,
        pssrb.grid_boundary_factors,
    )
end