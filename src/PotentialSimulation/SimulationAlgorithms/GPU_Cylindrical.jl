# function sor_cyl_gpu!(
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
    # # GPU Indexing
    # linear_idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x # CUDA.jl
    linear_idx = @index(Global) # KernelAbstractions.jl

    eff_size = broadcast(idim -> 2:size(potential, idim)-1, (1, 2, 3))
    if linear_idx <= prod(length.(eff_size))
        i1, i2, i3 = Tuple(CartesianIndices(eff_size)[linear_idx])
        # Comparison to CPU indecies: (Cyl / Car)
        # i3 <-> idx3 / ir / iz
        # i2 <-> idx2 / iφ / iy
        # i1 <-> idx1 / iz / ix
        in3 = i3 - 1
        in2 = i2 - 1
        in1 = nidx(i1, update_even_points, iseven(i2 + i3))
        ir = i3
        inr = in3
        iφ = i2
        inφ = in2
        iz = i1
        inz = in1
        izr = get_rbidx_right_neighbour(iz, update_even_points, iseven(i2 + i3))

        rb_tar_idx, rb_src_idx = update_even_points ? (rb_even::Int, rb_odd::Int) : (rb_odd::Int, rb_even::Int) 
                      
        pwwrr               = geom_weights[1][1, inr]
        pwwrl               = geom_weights[1][2, inr]
        r_inv_pwΔmpr        = geom_weights[1][3, inr]
        Δr_ext_inv_r_pwmprr = geom_weights[1][4, inr] 
        Δr_ext_inv_l_pwmprl = geom_weights[1][5, inr] 
        Δmpr_squared        = geom_weights[1][6, inr]  

        pwwφr        = geom_weights[2][1, inφ]
        pwwφl        = geom_weights[2][2, inφ]
        pwΔmpφ       = geom_weights[2][3, inφ]
        Δφ_ext_inv_r = geom_weights[2][4,  iφ]
        Δφ_ext_inv_l = geom_weights[2][4, inφ]

        if inr == 1
            pwwφr = T(0.5)
            pwwφl = T(0.5)
            pwΔmpφ = T(2π)
            Δφ_ext_inv_r = T(1/pwΔmpφ)
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

        pwwzr::T        = geom_weights[3][1, inz]
        pwwzl::T        = geom_weights[3][2, inz]
        pwΔmpz::T       = geom_weights[3][3, inz]
        Δz_ext_inv_r::T = geom_weights[3][4, inz + 1]
        Δz_ext_inv_l::T = geom_weights[3][4, inz]

        ϵ_rrr::T = ϵ_r[  ir,  iφ, inz + 1]
        ϵ_rlr::T = ϵ_r[  ir, inφ, inz + 1]
        ϵ_rrl::T = ϵ_r[  ir,  iφ, inz ]
        ϵ_rll::T = ϵ_r[  ir, inφ, inz ]
        ϵ_lrr::T = ϵ_r[ inr,  iφ, inz + 1]
        ϵ_llr::T = ϵ_r[ inr, inφ, inz + 1]
        ϵ_lrl::T = ϵ_r[ inr,  iφ, inz ] 
        ϵ_lll::T = ϵ_r[ inr, inφ, inz ] 

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

        old_potential::T = potential[iz, iφ, ir, rb_tar_idx]
        q_eff::T = is_weighting_potential ? zero(T) : (q_eff_imp[iz, iφ, ir, rb_tar_idx] + q_eff_fix[iz, iφ, ir, rb_tar_idx])

        wrr *= Δr_ext_inv_r_pwmprr_pwΔmpφ * pwΔmpz
        wrl *= Δr_ext_inv_l_pwmprl_pwΔmpφ * pwΔmpz
        wφr *= r_inv_pwΔmpr_Δφ_ext_inv_r * pwΔmpz
        wφl *= r_inv_pwΔmpr_Δφ_ext_inv_l * pwΔmpz
        wzr *= Δz_ext_inv_r * pwΔmpφ_Δmpr_squared
        wzl *= Δz_ext_inv_l * pwΔmpφ_Δmpr_squared

        vrr::T = potential[     iz,     iφ, ir + 1, rb_src_idx]
        vrl::T = potential[     iz,     iφ,    inr, rb_src_idx]
        vφr::T = only2d ? old_potential : potential[ iz, iφ + 1, ir, rb_src_idx]
        vφl::T = only2d ? old_potential : potential[ iz,    inφ, ir, rb_src_idx]
        vzr::T = potential[    izr,     iφ,     ir, rb_src_idx] 
        vzl::T = potential[izr - 1,     iφ,     ir, rb_src_idx]
        
        new_potential::T = calc_new_potential_SOR_3D(
            q_eff,
            volume_weights[iz, iφ, ir, rb_tar_idx],
            (wrr, wrl, wφr, wφl, wzr, wzl),
            (vrr, vrl, vφr, vφl, vzr, vzl),
            old_potential,
            sor_const[inr]
        )

        if depletion_handling_enabled
            new_potential, point_types[iz, iφ, ir, rb_tar_idx] = handle_depletion(
                new_potential,
                point_types[iz, iφ, ir, rb_tar_idx],
                (vrr, ifelse(inr == 1, vrr, vrl), vφr, vφl, vzr, vzl),
                q_eff_imp[iz, iφ, ir, rb_tar_idx],
                volume_weights[iz, iφ, ir, rb_tar_idx],
                sor_const[inr]
            )
        end

        potential[iz, iφ, ir, rb_tar_idx]::T = ifelse(point_types[iz, iφ, ir, rb_tar_idx] & update_bit > 0, new_potential, old_potential)
    end
end