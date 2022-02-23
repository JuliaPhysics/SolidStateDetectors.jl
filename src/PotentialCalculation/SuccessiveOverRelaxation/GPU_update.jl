@inline function update!(
    pssrb::PotentialCalculationSetup{T},
    kernel, 
    ndrange;
    use_nthreads::Int = Base.Threads.nthreads(),
    depletion_handling::Val{depletion_handling_enabled} = Val{false}(), 
    is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
    only2d::Val{only_2d} = Val{false}()
)::Nothing where {T, only_2d, depletion_handling_enabled, _is_weighting_potential}
    update_even_points = true
    wait(kernel( 
        pssrb.potential, pssrb.imp_scale, pssrb.point_types, pssrb.volume_weights, pssrb.q_eff_imp, pssrb.q_eff_fix, pssrb.ϵ_r,
        pssrb.geom_weights, pssrb.sor_const, update_even_points, depletion_handling_enabled, _is_weighting_potential, only_2d, 
        ndrange = ndrange
    ))
    apply_boundary_conditions!(pssrb, Val(update_even_points), only2d)
    update_even_points = false
    wait(kernel( 
        pssrb.potential, pssrb.imp_scale, pssrb.point_types, pssrb.volume_weights, pssrb.q_eff_imp, pssrb.q_eff_fix, pssrb.ϵ_r,
        pssrb.geom_weights, pssrb.sor_const, update_even_points, depletion_handling_enabled, _is_weighting_potential, only_2d,
        ndrange = ndrange
    ))
    apply_boundary_conditions!(pssrb, Val(update_even_points), only2d)
    return nothing
end
