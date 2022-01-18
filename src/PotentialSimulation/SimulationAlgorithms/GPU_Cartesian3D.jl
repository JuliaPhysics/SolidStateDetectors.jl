@kernel function sor_car_gpu!(
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
        iz = i3
        inz = in3
        iy = i2
        iny = in2
        ix = i1
        inx = in1
        ixr = get_rbidx_right_neighbour(ix, update_even_points, iseven(i2 + i3))
        
        rb_tar_idx, rb_src_idx = update_even_points ? (rb_even::Int, rb_odd::Int) : (rb_odd::Int, rb_even::Int) 

        pwwzr::T        = geom_weights[3][1, inz]
        pwwzl::T        = geom_weights[3][2, inz]
        pwΔmpz::T       = geom_weights[3][3, inz]
        Δz_ext_inv_r::T = geom_weights[3][4, inz + 1]
        Δz_ext_inv_l::T = geom_weights[3][4, inz]

        pwwyr::T  = geom_weights[2][1, iny]
        pwwyl::T  = geom_weights[2][2, iny]
        pwΔmpy::T = geom_weights[2][3, iny] 
        pwΔmpy_pwΔmpz::T = pwΔmpy * pwΔmpz
        Δy_ext_inv_r_pwΔmpz::T  = geom_weights[2][4, iny + 1] * pwΔmpz
        Δy_ext_inv_l_pwΔmpz::T  = geom_weights[2][4, iny]     * pwΔmpz
        Δz_ext_inv_r_pwΔmpy::T = Δz_ext_inv_r * pwΔmpy
        Δz_ext_inv_l_pwΔmpy::T = Δz_ext_inv_l * pwΔmpy

        pwwyr_pwwzr::T = pwwyr * pwwzr
        pwwyr_pwwzl::T = pwwyr * pwwzl
        pwwyl_pwwzr::T = pwwyl * pwwzr
        pwwyl_pwwzl::T = pwwyl * pwwzl

        pwwxr::T        = geom_weights[1][1, inx]
        pwwxl::T        = geom_weights[1][2, inx]
        pwΔmpx::T       = geom_weights[1][3, inx]
        Δx_ext_inv_r::T = geom_weights[1][4, inx + 1]
        Δx_ext_inv_l::T = geom_weights[1][4, inx]
        
        pwwxr_pwwzr::T = pwwxr * pwwzr
        pwwxl_pwwzr::T = pwwxl * pwwzr
        pwwxr_pwwzl::T = pwwxr * pwwzl
        pwwxl_pwwzl::T = pwwxl * pwwzl
        pwwxl_pwwyr::T = pwwxl * pwwyr
        pwwxr_pwwyr::T = pwwxr * pwwyr
        pwwxl_pwwyl::T = pwwxl * pwwyl
        pwwxr_pwwyl::T = pwwxr * pwwyl

        ϵ_rrr::T = ϵ_r[ inx + 1,  iy, iz]
        ϵ_rlr::T = ϵ_r[ inx + 1, iny, iz]
        ϵ_rrl::T = ϵ_r[ inx + 1,  iy, inz ]
        ϵ_rll::T = ϵ_r[ inx + 1, iny, inz ]
        ϵ_lrr::T = ϵ_r[ inx,  iy,  iz ]
        ϵ_llr::T = ϵ_r[ inx, iny,  iz ]
        ϵ_lrl::T = ϵ_r[ inx,  iy, inz ] 
        ϵ_lll::T = ϵ_r[ inx, iny, inz ] 

        # right weight in r: wrr
        wxr::T = ϵ_rrr * pwwyr_pwwzr
        wxr    = muladd(ϵ_rlr, pwwyl_pwwzr, wxr)   
        wxr    = muladd(ϵ_rrl, pwwyr_pwwzl, wxr)    
        wxr    = muladd(ϵ_rll, pwwyl_pwwzl, wxr)
        # left weight in r: wrr
        wxl::T = ϵ_lrr * pwwyr_pwwzr
        wxl    = muladd(ϵ_llr, pwwyl_pwwzr, wxl)   
        wxl    = muladd(ϵ_lrl, pwwyr_pwwzl, wxl)    
        wxl    = muladd(ϵ_lll, pwwyl_pwwzl, wxl) 
        # right weight in φ: wφr
        wyr::T = ϵ_lrr * pwwxl_pwwzr 
        wyr    = muladd(ϵ_rrr, pwwxr_pwwzr, wyr)  
        wyr    = muladd(ϵ_lrl, pwwxl_pwwzl, wyr)    
        wyr    = muladd(ϵ_rrl, pwwxr_pwwzl, wyr) 
        # left weight in φ: wφl
        wyl::T = ϵ_llr * pwwxl_pwwzr 
        wyl    = muladd(ϵ_rlr, pwwxr_pwwzr, wyl)  
        wyl    = muladd(ϵ_lll, pwwxl_pwwzl, wyl)    
        wyl    = muladd(ϵ_rll, pwwxr_pwwzl, wyl) 
        # right weight in z: wzr
        wzr::T = ϵ_rrr * pwwxr_pwwyr  
        wzr    = muladd(ϵ_rlr, pwwxr_pwwyl, wzr)     
        wzr    = muladd(ϵ_lrr, pwwxl_pwwyr, wzr)     
        wzr    = muladd(ϵ_llr, pwwxl_pwwyl, wzr)
        # left weight in z: wzr
        wzl::T = ϵ_rrl * pwwxr_pwwyr 
        wzl    = muladd(ϵ_rll, pwwxr_pwwyl, wzl)    
        wzl    = muladd(ϵ_lrl, pwwxl_pwwyr, wzl)    
        wzl    = muladd(ϵ_lll, pwwxl_pwwyl, wzl)

        old_potential::T = potential[ix, iy, iz, rb_tar_idx]
        q_eff::T = is_weighting_potential ? zero(T) : (q_eff_imp[ix, iy, iz, rb_tar_idx] + q_eff_fix[ix, iy, iz, rb_tar_idx])
        
        wxr *= pwΔmpy_pwΔmpz * Δx_ext_inv_r 
        wxl *= pwΔmpy_pwΔmpz * Δx_ext_inv_l
        wyr *= Δy_ext_inv_r_pwΔmpz * pwΔmpx
        wyl *= Δy_ext_inv_l_pwΔmpz * pwΔmpx
        wzr *= Δz_ext_inv_r_pwΔmpy * pwΔmpx 
        wzl *= Δz_ext_inv_l_pwΔmpy * pwΔmpx 

        vxr::T = potential[    ixr,     iy,     iz, rb_src_idx] 
        vxl::T = potential[ixr - 1,     iy,     iz, rb_src_idx]
        vyr::T = potential[     ix, iy + 1,     iz, rb_src_idx]
        vyl::T = potential[     ix,    iny,     iz, rb_src_idx]
        vzr::T = potential[     ix,     iy, iz + 1, rb_src_idx]
        vzl::T = potential[     ix,     iy,    inz, rb_src_idx]

        new_potential::T = calc_new_potential_SOR_3D(
            q_eff,
            volume_weights[ix, iy, iz, rb_tar_idx],
            (wxr, wxl, wyr, wyl, wzr, wzl),
            (vxr, vxl, vyr, vyl, vzr, vzl),
            old_potential,
            sor_const[1]
        )

        if depletion_handling_enabled
            vmin::T = min(vxr, vxl, vyr, vyl, vzr, vzl)
            vmax::T = max(vxr, vxl, vyr, vyl, vzr, vzl)
            
            new_point_type = point_types[ix, iy, iz, rb_tar_idx] 
            if new_potential < vmin || new_potential > vmax
                new_potential -= q_eff_imp[ix, iy, iz, rb_tar_idx] * volume_weights[ix, iy, iz, rb_tar_idx] * sor_const[1]
                if (point_types[ix, iy, iz, rb_tar_idx] & undepleted_bit == 0 &&
                    point_types[ix, iy, iz, rb_tar_idx] & pn_junction_bit > 0) 
                    new_point_type += undepleted_bit
                end # mark this point as undepleted
            elseif point_types[ix, iy, iz, rb_tar_idx] & undepleted_bit > 0
                new_point_type -= undepleted_bit
            end
            point_types[ix, iy, iz, rb_tar_idx] = new_point_type
        end 

        potential[ix, iy, iz, rb_tar_idx]::T = ifelse(point_types[ix, iy, iz, rb_tar_idx] & update_bit > 0, new_potential, old_potential)
    end
end