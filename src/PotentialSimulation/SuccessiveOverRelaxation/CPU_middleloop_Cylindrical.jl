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

function load_weights_for_innerloop!(line_weights, pssrb::PotentialSimulationSetupRB{T, Cylindrical, 3, Array{T, 3}},
    iφ, inφ, ir, inr,
    update_even_points, rφi_is_even_t::Val{rφi_is_even},
    pwwrr, pwwrl, pwwφr, pwwφl, 
    pwwrr_pwwφr, pwwrl_pwwφr, pwwrr_pwwφl, pwwrl_pwwφl,
    pwΔmpφ_Δmpr_squared, 
    Δr_ext_inv_r_pwmprr_pwΔmpφ, Δr_ext_inv_l_pwmprl_pwΔmpφ, 
    r_inv_pwΔmpr_Δφ_ext_inv_r, r_inv_pwΔmpr_Δφ_ext_inv_l,
    ) where {T, rφi_is_even}
# Prepare weights for most inner loop of the SOR
@fastmath @inbounds @simd ivdep for iz in 2:(size(pssrb.potential, 1) - 1)
inz::Int = nidx(iz, update_even_points, rφi_is_even_t)

pwwzr::T        = pssrb.geom_weights[3][1, inz]
pwwzl::T        = pssrb.geom_weights[3][2, inz]
pwΔmpz::T       = pssrb.geom_weights[3][3, inz]
Δz_ext_inv_r::T = pssrb.geom_weights[3][4, inz + 1]
Δz_ext_inv_l::T = pssrb.geom_weights[3][4, inz]

ϵ_rrr::T = pssrb.ϵ_r[  ir,  iφ, inz + 1]
ϵ_rlr::T = pssrb.ϵ_r[  ir, inφ, inz + 1]
ϵ_rrl::T = pssrb.ϵ_r[  ir,  iφ, inz ]
ϵ_rll::T = pssrb.ϵ_r[  ir, inφ, inz ]
ϵ_lrr::T = pssrb.ϵ_r[ inr,  iφ, inz + 1]
ϵ_llr::T = pssrb.ϵ_r[ inr, inφ, inz + 1]
ϵ_lrl::T = pssrb.ϵ_r[ inr,  iφ, inz ] 
ϵ_lll::T = pssrb.ϵ_r[ inr, inφ, inz ] 

pwwφr_pwwzr::T = pwwφr * pwwzr
pwwφl_pwwzr::T = pwwφl * pwwzr
pwwφr_pwwzl::T = pwwφr * pwwzl
pwwφl_pwwzl::T = pwwφl * pwwzl
pwwrl_pwwzr::T = pwwrl * pwwzr
pwwrr_pwwzr::T = pwwrr * pwwzr
pwwrl_pwwzl::T = pwwrl * pwwzl
pwwrr_pwwzl::T = pwwrr * pwwzl

# right weight in r: wrr
wrr::T = ϵ_rrr * pwwφr_pwwzr
wrr    = muladd(ϵ_rlr, pwwφl_pwwzr, wrr)   
wrr    = muladd(ϵ_rrl, pwwφr_pwwzl, wrr)    
wrr    = muladd(ϵ_rll, pwwφl_pwwzl, wrr)
# left weight in r: wrr
wrl::T = ϵ_lrr * pwwφr_pwwzr
wrl    = muladd(ϵ_llr, pwwφl_pwwzr, wrl)   
wrl    = muladd(ϵ_lrl, pwwφr_pwwzl, wrl)    
wrl    = muladd(ϵ_lll, pwwφl_pwwzl, wrl) 
# right weight in φ: wφr
wφr::T = ϵ_lrr * pwwrl_pwwzr 
wφr    = muladd(ϵ_rrr, pwwrr_pwwzr, wφr)  
wφr    = muladd(ϵ_lrl, pwwrl_pwwzl, wφr)    
wφr    = muladd(ϵ_rrl, pwwrr_pwwzl, wφr) 
# left weight in φ: wφl
wφl::T = ϵ_llr * pwwrl_pwwzr 
wφl    = muladd(ϵ_rlr, pwwrr_pwwzr, wφl)  
wφl    = muladd(ϵ_lll, pwwrl_pwwzl, wφl)    
wφl    = muladd(ϵ_rll, pwwrr_pwwzl, wφl) 
# right weight in z: wzr
wzr::T = ϵ_rrr * pwwrr_pwwφr  
wzr    = muladd(ϵ_rlr, pwwrr_pwwφl, wzr)     
wzr    = muladd(ϵ_lrr, pwwrl_pwwφr, wzr)     
wzr    = muladd(ϵ_llr, pwwrl_pwwφl, wzr)
# left weight in z: wzr
wzl::T = ϵ_rrl * pwwrr_pwwφr 
wzl    = muladd(ϵ_rll, pwwrr_pwwφl, wzl)    
wzl    = muladd(ϵ_lrl, pwwrl_pwwφr, wzl)    
wzl    = muladd(ϵ_lll, pwwrl_pwwφl, wzl)

wrl *= Δr_ext_inv_l_pwmprl_pwΔmpφ * pwΔmpz
wrr *= Δr_ext_inv_r_pwmprr_pwΔmpφ * pwΔmpz
wφl *= r_inv_pwΔmpr_Δφ_ext_inv_l * pwΔmpz
wφr *= r_inv_pwΔmpr_Δφ_ext_inv_r * pwΔmpz
wzl *= Δz_ext_inv_l * pwΔmpφ_Δmpr_squared
wzr *= Δz_ext_inv_r * pwΔmpφ_Δmpr_squared

line_weights[iz-1, 1] = wrl 
line_weights[iz-1, 2] = wrr 
line_weights[iz-1, 3] = wφl 
line_weights[iz-1, 4] = wφr 
line_weights[iz-1, 5] = wzl 
line_weights[iz-1, 6] = wzr 
end # z loop
nothing
end