@inline function get_geom_weights_outerloop(
    geom_weights::NTuple{3, <:AbstractArray{T, 2}}, i, ::Type{Cylindrical}
) where {T}
    geom_weights[1][1, i],
    geom_weights[1][2, i],
    geom_weights[1][3, i],
    geom_weights[1][4, i], 
    geom_weights[1][5, i], 
    geom_weights[1][6, i]  
end

@inline function get_geom_weights_outerloop(
    geom_weights::NTuple{3, <:AbstractArray{T, 2}}, i, ::Type{Cartesian}
) where {T}
    geom_weights[1][1, i],
    geom_weights[1][2, i],
    geom_weights[1][3, i],
    geom_weights[1][4, i + 1],
    geom_weights[1][4, i]
end

@inline function prepare_weights_in_middleloop(
    geom_weights::NTuple{3, <:AbstractArray{T, 2}}, ::Type{Cylindrical},
    i2, in2,
    pwwrr, pwwrl, r_inv_pwΔmpr, Δr_ext_inv_r_pwmprr, Δr_ext_inv_l_pwmprl, Δmpr_squared, 
    is_r0_t::Val{is_r0}
) where {T, is_r0}
    pwwφr        = geom_weights[2][1, in2]
    pwwφl        = geom_weights[2][2, in2]
    pwΔmpφ       = geom_weights[2][3, in2]
    Δφ_ext_inv_r = geom_weights[2][4,  i2]
    Δφ_ext_inv_l = geom_weights[2][4, in2]

    if is_r0
        pwwφr = T(0.5)
        pwwφl = T(0.5)
        pwΔmpφ = T(2π)
        Δφ_ext_inv_r = inv(pwΔmpφ)
        Δφ_ext_inv_l = Δφ_ext_inv_r
    end
    pwwrr_pwwφr = pwwrr * pwwφr
    pwwrl_pwwφr = pwwrl * pwwφr
    pwwrr_pwwφl = pwwrr * pwwφl
    pwwrl_pwwφl = pwwrl * pwwφl

    pwΔmpφ_Δmpr_squared = pwΔmpφ * Δmpr_squared
    Δr_ext_inv_r_pwmprr_pwΔmpφ = Δr_ext_inv_r_pwmprr * pwΔmpφ
    Δr_ext_inv_l_pwmprl_pwΔmpφ = Δr_ext_inv_l_pwmprl * pwΔmpφ
    r_inv_pwΔmpr_Δφ_ext_inv_r = r_inv_pwΔmpr * Δφ_ext_inv_r
    r_inv_pwΔmpr_Δφ_ext_inv_l = r_inv_pwΔmpr * Δφ_ext_inv_l
    return (
        pwwrr, pwwrl, pwwφr, pwwφl, 
        pwwrr_pwwφr, pwwrl_pwwφr, pwwrr_pwwφl, pwwrl_pwwφl,
        pwΔmpφ_Δmpr_squared,
        Δr_ext_inv_r_pwmprr_pwΔmpφ, Δr_ext_inv_l_pwmprl_pwΔmpφ,
        r_inv_pwΔmpr_Δφ_ext_inv_r, r_inv_pwΔmpr_Δφ_ext_inv_l
    )
end

@inline function prepare_weights_in_middleloop(
    geom_weights::NTuple{3, <:AbstractArray{T, 2}}, ::Type{Cartesian},
    i2, in2,
    pww3r, pww3l, pwΔmp3, Δ3_ext_inv_r, Δ3_ext_inv_l,
    is_r0_t
) where {T}
    pww2r  = geom_weights[2][1, in2]
    pww2l  = geom_weights[2][2, in2]
    pwΔmp2 = geom_weights[2][3, in2] 
    pwΔmp2_pwΔmp3 = pwΔmp2 * pwΔmp3
    Δ2_ext_inv_r_pwΔmp3  = geom_weights[2][4, in2 + 1] * pwΔmp3
    Δ2_ext_inv_l_pwΔmp3  = geom_weights[2][4, in2]     * pwΔmp3
    Δ3_ext_inv_r_pwΔmp2 = Δ3_ext_inv_r * pwΔmp2
    Δ3_ext_inv_l_pwΔmp2 = Δ3_ext_inv_l * pwΔmp2

    pww2r_pww3r = pww2r * pww3r
    pww2r_pww3l = pww2r * pww3l
    pww2l_pww3r = pww2l * pww3r
    pww2l_pww3l = pww2l * pww3l
    return (
        pww3r, pww3l, pww2r, pww2l,
        pww2r_pww3r, pww2r_pww3l, pww2l_pww3r, pww2l_pww3l,
        pwΔmp2_pwΔmp3,
        Δ3_ext_inv_r_pwΔmp2, Δ3_ext_inv_l_pwΔmp2,
        Δ2_ext_inv_r_pwΔmp3, Δ2_ext_inv_l_pwΔmp3, 
    )
end

@fastmath function middleloop!( 
    i3::Int, rb_tar_idx::Int, rb_src_idx::Int, 
    pcs::PotentialCalculationSetup{T, S},
    update_even_points::Val{even_points},
    depletion_handling::Val{depletion_handling_enabled},
    is_weighting_potential::Val{_is_weighting_potential}, 
    only2d::Val{only_2d},
    idx3iseven::Val{idx3_is_even},
    is_r0_t::Val{is_r0}
)::Nothing where {T, S, even_points, depletion_handling_enabled, _is_weighting_potential, only_2d, idx3_is_even, is_r0}
    @inbounds begin 
        in3 = i3 - 1 
                
        geom_weights_3 = get_geom_weights_outerloop(pcs.geom_weights, in3, S)

        line_weights::Array{T, 2} = Array{T, 2}(undef, size(pcs.potential, 1) - 2, 6)
        # Even though this causes some allocations it 
        # is faster than using a predefined array, e.g. stored in pcs
        # Especially when using multiple threads

        #=
            Splitting this into two loops over even and uneven i2
            seems to increase the performance quite a bit by hard-coding oscillating type 
            instead of one loop using:
                rφi_is_even_t::Union{Val{true}, Val{false}} = Val(iseven(i3 + i2) ? true : false)
        =#
        for i2 in (only_2d ? (2,) : 2:2:(size(pcs.potential, 2) - 1))
            in2 = i2 - 1
            i23_is_even_t = Val(idx3_is_even ? true : false)

            geom_weights_2 = prepare_weights_in_middleloop(
                pcs.geom_weights, S, i2, in2, 
                geom_weights_3...,
                is_r0_t
            )

            calculate_weights_for_innerloop!(line_weights, pcs, i2, in2, i3, in3,
                update_even_points, i23_is_even_t, 
                geom_weights_2...
            )

            innerloop!(line_weights, pcs, i2, in2, i3, in3, rb_tar_idx, rb_src_idx, 
                update_even_points, i23_is_even_t, 
                depletion_handling, is_weighting_potential, only2d)
        end 
        for i2 in 3:2:(size(pcs.potential, 2) - 1)
            in2 = i2 - 1
            i23_is_even_t = Val(idx3_is_even ? false : true)

            geom_weights_2 = prepare_weights_in_middleloop(
                pcs.geom_weights, S, i2, in2, 
                geom_weights_3...,
                is_r0_t
            )
            
            calculate_weights_for_innerloop!(line_weights, pcs, i2, in2, i3, in3,
                update_even_points, i23_is_even_t, 
                geom_weights_2...
            )

            innerloop!(line_weights, pcs, i2, in2, i3, in3, rb_tar_idx, rb_src_idx, 
                update_even_points, i23_is_even_t, 
                depletion_handling, is_weighting_potential, only2d)
        end
    end
end

