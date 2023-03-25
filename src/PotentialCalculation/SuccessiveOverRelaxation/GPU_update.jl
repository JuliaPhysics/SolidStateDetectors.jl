@inline function update!(
    pcs::PotentialCalculationSetup{T},
    kernel, 
    ndrange;
    use_nthreads::Int = Base.Threads.nthreads(),
    depletion_handling::Val{depletion_handling_enabled} = Val{false}(), 
    is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
    only2d::Val{only_2d} = Val{false}()
)::Nothing where {T, only_2d, depletion_handling_enabled, _is_weighting_potential}
    update_even_points = true
    _ka_synchronize(kernel, kernel( 
        pcs.potential, pcs.imp_scale, pcs.point_types, pcs.volume_weights, pcs.q_eff_imp, pcs.q_eff_fix, pcs.ϵ_r,
        pcs.geom_weights, pcs.sor_const, update_even_points, depletion_handling_enabled, _is_weighting_potential, only_2d, 
        ndrange = ndrange
    ))
    apply_boundary_conditions!(pcs, Val(update_even_points), only2d)
    update_even_points = false
    _ka_synchronize(kernel, kernel( 
        pcs.potential, pcs.imp_scale, pcs.point_types, pcs.volume_weights, pcs.q_eff_imp, pcs.q_eff_fix, pcs.ϵ_r,
        pcs.geom_weights, pcs.sor_const, update_even_points, depletion_handling_enabled, _is_weighting_potential, only_2d,
        ndrange = ndrange
    ))
    apply_boundary_conditions!(pcs, Val(update_even_points), only2d)
    return nothing
end
