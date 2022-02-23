@kernel function sor_cyl_gpu!(
    potential::AbstractArray{T, 4},
    imp_scale::AbstractArray{T, 4},
    @Const(point_types::AbstractArray{PointType, 4}),
    @Const(volume_weights::AbstractArray{T, 4}),
    @Const(q_eff_imp::AbstractArray{T, 4}),
    @Const(q_eff_fix::AbstractArray{T, 4}),
    @Const(ϵ_r::AbstractArray{T, 3}),
    @Const(geom_weights::NTuple{3, AbstractArray{T, 2}}),
    @Const(sor_const::AbstractArray{T, 1}),
    @Const(update_even_points::Bool),
    @Const(depletion_handling_enabled::Bool),
    @Const(is_weighting_potential::Bool),
    @Const(only2d::Bool)
) where {T}
    linear_idx = @index(Global)     
    sor_kernel(
        potential,
        imp_scale,
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
