@kernel function sor_cyl_gpu!(
    potential::AbstractArray{T, 4},
    point_types::AbstractArray{PointType, 4},
    volume_weights::AbstractArray{T, 4},
    q_eff_imp::AbstractArray{T, 4},
    q_eff_fix::AbstractArray{T, 4},
    ϵ_r::AbstractArray{T, 3},
    geom_weights::NTuple{3, <:AbstractArray{T, 2}},
    sor_const::AbstractArray{T, 1},
    update_even_points::Bool,
    depletion_handling_enabled::Bool,
    is_weighting_potential::Bool,
    only2d::Bool
) where {T}
    linear_idx = @index(Global)     
    sor_kernel(
        potential,
        point_types,
        volume_weights,
        q_eff_imp,
        q_eff_fix,
        ϵ_r,
        geom_weights,
        sor_const,
        update_even_points,
        depletion_handling_enabled,
        is_weighting_potential,
        only2d, 
        Cylindrical,
        linear_idx
    )
end
