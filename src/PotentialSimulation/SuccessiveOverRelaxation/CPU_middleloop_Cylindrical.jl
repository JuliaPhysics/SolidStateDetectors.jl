@fastmath function middleloop!( ir::Int, rb_tar_idx::Int, rb_src_idx::Int, 
                                pssrb::PotentialSimulationSetupRB{T, Cylindrical, 3, Array{T, 3}},
                                update_even_points::Val{even_points},
                                depletion_handling::Val{depletion_handling_enabled},
                                is_weighting_potential::Val{_is_weighting_potential}, 
                                only2d::Val{only_2d},
                                idx3iseven::Val{idx3_is_even})::Nothing where {T, even_points, depletion_handling_enabled, _is_weighting_potential, only_2d, idx3_is_even}
    @inbounds begin 
        inr = ir - 1 
                
        pwwrr               = pssrb.geom_weights[1][1, inr]
        pwwrl               = pssrb.geom_weights[1][2, inr]
        r_inv_pwΔmpr        = pssrb.geom_weights[1][3, inr]
        Δr_ext_inv_r_pwmprr = pssrb.geom_weights[1][4, inr] 
        Δr_ext_inv_l_pwmprl = pssrb.geom_weights[1][5, inr] 
        Δmpr_squared        = pssrb.geom_weights[1][6, inr]  

        line_weights::Array{T, 2} = Array{T, 2}(undef, size(pssrb.potential, 1) - 2, 6)
        # Even though this causes some allocations it is faster than using a predefined array, e.g. stored in pssrb
        # Especially when using multiple threads

        if only_2d
            iφ = 2
            inφ = iφ - 1

            pwwφr        = pssrb.geom_weights[2][1, inφ]
            pwwφl        = pssrb.geom_weights[2][2, inφ]
            pwΔmpφ       = pssrb.geom_weights[2][3, inφ]
            Δφ_ext_inv_r = pssrb.geom_weights[2][4,  iφ]
            Δφ_ext_inv_l = pssrb.geom_weights[2][4, inφ]

            if inr == 1
                pwwφr = T(0.5)
                pwwφl = T(0.5)
                pwΔmpφ = T(2π)
                Δφ_ext_inv_r = inv(pwΔmpφ)
                Δφ_ext_inv_l = Δφ_ext_inv_r
            end
            pwwrr_pwwφr = pwwrr * pwwφr
            pwwrr_pwwφl = pwwrr * pwwφl
            pwwrl_pwwφr = pwwrl * pwwφr
            pwwrl_pwwφl = pwwrl * pwwφl

            Δr_ext_inv_r_pwmprr_pwΔmpφ = Δr_ext_inv_r_pwmprr * pwΔmpφ
            Δr_ext_inv_l_pwmprl_pwΔmpφ = Δr_ext_inv_l_pwmprl * pwΔmpφ
            pwΔmpφ_Δmpr_squared = pwΔmpφ * Δmpr_squared
            r_inv_pwΔmpr_Δφ_ext_inv_r = r_inv_pwΔmpr * Δφ_ext_inv_r
            r_inv_pwΔmpr_Δφ_ext_inv_l = r_inv_pwΔmpr * Δφ_ext_inv_l

            load_weights_for_innerloop!(line_weights, pssrb, iφ, inφ, ir, inr,
                update_even_points, Val(idx3_is_even ? true : false), 
                pwwrr, pwwrl, pwwφr, pwwφl, 
                pwwrr_pwwφr, pwwrl_pwwφr, pwwrr_pwwφl, pwwrl_pwwφl,
                pwΔmpφ_Δmpr_squared, 
                Δr_ext_inv_r_pwmprr_pwΔmpφ, Δr_ext_inv_l_pwmprl_pwΔmpφ, 
                r_inv_pwΔmpr_Δφ_ext_inv_r, r_inv_pwΔmpr_Δφ_ext_inv_l,
            )

            innerloop!(line_weights, pssrb, iφ, inφ, ir, inr, rb_tar_idx, rb_src_idx, 
                update_even_points, Val(idx3_is_even ? true : false), 
                depletion_handling, is_weighting_potential, only2d)
        else
            #=
                Splitting this into two loops over even and uneven iφ
                seems to increase the performance quite a bit
                by hard-coding oscillating type: 
                    rφi_is_even_t::Union{Val{true}, Val{false}} = Val(iseven(ir + iφ) ? true : false)
            =#
            for iφ in 2:2:(size(pssrb.potential, 2) - 1)
                inφ = iφ - 1

                pwwφr        = pssrb.geom_weights[2][1, inφ]
                pwwφl        = pssrb.geom_weights[2][2, inφ]
                pwΔmpφ       = pssrb.geom_weights[2][3, inφ]
                Δφ_ext_inv_r = pssrb.geom_weights[2][4,  iφ]
                Δφ_ext_inv_l = pssrb.geom_weights[2][4, inφ]

                if inr == 1
                    pwwφr = T(0.5)
                    pwwφl = T(0.5)
                    pwΔmpφ = T(2π)
                    Δφ_ext_inv_r = inv(pwΔmpφ)
                    Δφ_ext_inv_l = Δφ_ext_inv_r
                end
                pwwrr_pwwφr = pwwrr * pwwφr
                pwwrr_pwwφl = pwwrr * pwwφl
                pwwrl_pwwφr = pwwrl * pwwφr
                pwwrl_pwwφl = pwwrl * pwwφl

                Δr_ext_inv_r_pwmprr_pwΔmpφ = Δr_ext_inv_r_pwmprr * pwΔmpφ
                Δr_ext_inv_l_pwmprl_pwΔmpφ = Δr_ext_inv_l_pwmprl * pwΔmpφ
                pwΔmpφ_Δmpr_squared = pwΔmpφ * Δmpr_squared
                r_inv_pwΔmpr_Δφ_ext_inv_r = r_inv_pwΔmpr * Δφ_ext_inv_r
                r_inv_pwΔmpr_Δφ_ext_inv_l = r_inv_pwΔmpr * Δφ_ext_inv_l

                load_weights_for_innerloop!(line_weights, pssrb, iφ, inφ, ir, inr,
                    update_even_points, Val(idx3_is_even ? true : false), 
                    pwwrr, pwwrl, pwwφr, pwwφl, 
                    pwwrr_pwwφr, pwwrl_pwwφr, pwwrr_pwwφl, pwwrl_pwwφl,
                    pwΔmpφ_Δmpr_squared, 
                    Δr_ext_inv_r_pwmprr_pwΔmpφ, Δr_ext_inv_l_pwmprl_pwΔmpφ, 
                    r_inv_pwΔmpr_Δφ_ext_inv_r, r_inv_pwΔmpr_Δφ_ext_inv_l,
                )

                innerloop!(line_weights, pssrb, iφ, inφ, ir, inr, rb_tar_idx, rb_src_idx, 
                    update_even_points, Val(idx3_is_even ? true : false), 
                    depletion_handling, is_weighting_potential, only2d)
            end 
            for iφ in 3:2:(size(pssrb.potential, 2) - 1)
                inφ = iφ - 1
                
                pwwφr        = pssrb.geom_weights[2][1, inφ]
                pwwφl        = pssrb.geom_weights[2][2, inφ]
                pwΔmpφ       = pssrb.geom_weights[2][3, inφ]
                Δφ_ext_inv_r = pssrb.geom_weights[2][4,  iφ]
                Δφ_ext_inv_l = pssrb.geom_weights[2][4, inφ]

                if inr == 1
                    pwwφr = T(0.5)
                    pwwφl = T(0.5)
                    pwΔmpφ = T(2π)
                    Δφ_ext_inv_r = inv(pwΔmpφ)
                    Δφ_ext_inv_l = Δφ_ext_inv_r
                end
                pwwrr_pwwφr = pwwrr * pwwφr
                pwwrr_pwwφl = pwwrr * pwwφl
                pwwrl_pwwφr = pwwrl * pwwφr
                pwwrl_pwwφl = pwwrl * pwwφl

                Δr_ext_inv_r_pwmprr_pwΔmpφ = Δr_ext_inv_r_pwmprr * pwΔmpφ
                Δr_ext_inv_l_pwmprl_pwΔmpφ = Δr_ext_inv_l_pwmprl * pwΔmpφ
                pwΔmpφ_Δmpr_squared = pwΔmpφ * Δmpr_squared
                r_inv_pwΔmpr_Δφ_ext_inv_r = r_inv_pwΔmpr * Δφ_ext_inv_r
                r_inv_pwΔmpr_Δφ_ext_inv_l = r_inv_pwΔmpr * Δφ_ext_inv_l

                load_weights_for_innerloop!(line_weights, pssrb, iφ, inφ, ir, inr,
                    update_even_points, Val(idx3_is_even ? false : true), 
                    pwwrr, pwwrl, pwwφr, pwwφl,  
                    pwwrr_pwwφr, pwwrl_pwwφr, pwwrr_pwwφl, pwwrl_pwwφl,
                    pwΔmpφ_Δmpr_squared, 
                    Δr_ext_inv_r_pwmprr_pwΔmpφ, Δr_ext_inv_l_pwmprl_pwΔmpφ, 
                    r_inv_pwΔmpr_Δφ_ext_inv_r, r_inv_pwΔmpr_Δφ_ext_inv_l,
                )

                innerloop!(line_weights, pssrb, iφ, inφ, ir, inr, rb_tar_idx, rb_src_idx, 
                    update_even_points, Val(idx3_is_even ? false : true), 
                    depletion_handling, is_weighting_potential, only2d)
            end # φ loop            
        end
    end # inbounds
end

