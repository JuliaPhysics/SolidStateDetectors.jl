function innerloop!(
    line_weights, pssrb::PotentialSimulationSetupRB{T, S, 3, Array{T, 3}}, 
    i2, in2, i3, in3, rb_tar_idx, rb_src_idx, 
    update_even_points, i23_is_even_t,
    depletion_handling::Val{depletion_handling_enabled},
    is_weighting_potential::Val{_is_weighting_potential}, 
    only2d::Val{only_2d}
) where {T, S, depletion_handling_enabled, _is_weighting_potential, only_2d}
    @fastmath @inbounds @simd ivdep for i1 in 2:(size(pssrb.potential, 1) - 1)
        i1r = get_rbidx_right_neighbour(i1, update_even_points, i23_is_even_t)
        
        old_potential = pssrb.potential[i1, i2, i3, rb_tar_idx]
        q_eff = _is_weighting_potential ? zero(T) : (pssrb.q_eff_imp[i1, i2, i3, rb_tar_idx] + pssrb.q_eff_fix[i1, i2, i3, rb_tar_idx])

        weights = get_sor_weights(line_weights, i1-1)

        neighbor_potentials = get_neighbor_potentials(
            pssrb, old_potential, i1, i2, i3, i1r, in2, in3, rb_src_idx, only2d
        )
        
        new_potential = calc_new_potential_SOR_3D(
            q_eff,
            pssrb.volume_weights[i1, i2, i3, rb_tar_idx],
            weights,
            neighbor_potentials,
            old_potential,
            get_sor_constant(pssrb, in3)
        )

        if depletion_handling_enabled
            new_potential, pssrb.point_types[i1, i2, i3, rb_tar_idx] = handle_depletion(
                new_potential,
                pssrb.point_types[i1, i2, i3, rb_tar_idx],
                r0_handling_depletion_handling(neighbor_potentials, S, in3),
                pssrb.q_eff_imp[i1, i2, i3, rb_tar_idx],
                pssrb.volume_weights[i1, i2, i3, rb_tar_idx],
                get_sor_constant(pssrb, in3)
            )
        end

        pssrb.potential[i1, i2, i3, rb_tar_idx] = ifelse(pssrb.point_types[i1, i2, i3, rb_tar_idx] & update_bit > 0, new_potential, old_potential)
    end
    nothing
end

@inline function get_sor_weights(line_weights::Matrix{T}, i::Int)::NTuple{6, T} where {T}
    @inbounds return (
        line_weights[i, 1], 
        line_weights[i, 2], 
        line_weights[i, 3], 
        line_weights[i, 4], 
        line_weights[i, 5], 
        line_weights[i, 6]
    ) 
end

@inline function get_neighbor_potentials(
    pssrb::PotentialSimulationSetupRB{T, S, 3, Array{T, 3}},
    old_potential, i1, i2, i3, i1r, in2, in3, rb_src_idx, only2d::Val{only_2d}
)::NTuple{6, T} where {T,S,only_2d}
    @inbounds return (
        pssrb.potential[     i1,     i2,    in3, rb_src_idx],
        pssrb.potential[     i1,     i2, i3 + 1, rb_src_idx],
        only_2d ? old_potential : pssrb.potential[ i1,    in2, i3, rb_src_idx],
        only_2d ? old_potential : pssrb.potential[ i1, i2 + 1, i3, rb_src_idx],
        pssrb.potential[i1r - 1,     i2,     i3, rb_src_idx],
        pssrb.potential[    i1r,     i2,     i3, rb_src_idx]
    ) 
end

@inline function get_sor_constant(pssrb::PotentialSimulationSetupRB{T, Cylindrical}, i::Int)::T where {T}
    @inbounds pssrb.sor_const[i]
end
@inline function get_sor_constant(pssrb::PotentialSimulationSetupRB{T, Cartesian}, i::Int)::T where {T}
    @inbounds pssrb.sor_const[1]
end

@inline function r0_handling_depletion_handling(
    np::NTuple{6, T}, ::Type{Cylindrical}, i::Int
) where {T}
    return (ifelse(i == 1, np[2], np[1]), np[2:6]...)
end
@inline function r0_handling_depletion_handling(
    np::NTuple{6, T}, ::Type{Cartesian}, i::Int
) where {T}
    return np
end