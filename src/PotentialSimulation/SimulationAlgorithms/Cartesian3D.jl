@fastmath function middleloop!( iz::Int, rb_tar_idx::Int, rb_src_idx::Int, 
                                pssrb::PotentialSimulationSetupRB{T, Cartesian, 3, Array{T, 3}},
                                update_even_points::Val{even_points},
                                depletion_handling::Val{depletion_handling_enabled},
                                is_weighting_potential::Val{_is_weighting_potential},
                                only2d::Val{only_2d},
                                idx3iseven::Val{idx3_is_even})::Nothing where {T, even_points, depletion_handling_enabled, _is_weighting_potential, only_2d, idx3_is_even}
    @inbounds begin 
        inz::Int = iz - 1 
                
        pwwzr::T        = pssrb.geom_weights[3][1, inz]
        pwwzl::T        = pssrb.geom_weights[3][2, inz]
        pwΔmpz::T       = pssrb.geom_weights[3][3, inz]
        Δz_ext_inv_r::T = pssrb.geom_weights[3][4, inz + 1]
        Δz_ext_inv_l::T = pssrb.geom_weights[3][4, inz]

        line_weights::Array{T, 2} = Array{T, 2}(undef, size(pssrb.potential, 1) - 2, 6)
        # Even though this causes some allocations it is faster than using a predefined array, e.g. stored in pssrb
        # Especially when using multiple threads
        
        #=
            Splitting this into two loops over even and uneven iφ
            seems to increase the performance quite a bit
            by hard-coding oscillating type: 
                zyi_is_even_t::Union{Val{true}, Val{false}} = Val(iseven(iz + iy) ? true : false)
        =#
        for iy in 2:2:(size(pssrb.potential, 2) - 1)
            iny::Int = iy - 1

            pwwyr::T  = pssrb.geom_weights[2][1, iny]
            pwwyl::T  = pssrb.geom_weights[2][2, iny]
            pwΔmpy::T = pssrb.geom_weights[2][3, iny] 
            pwΔmpy_pwΔmpz::T = pwΔmpy * pwΔmpz
            Δy_ext_inv_r_pwΔmpz::T  = pssrb.geom_weights[2][4, iny + 1] * pwΔmpz
            Δy_ext_inv_l_pwΔmpz::T  = pssrb.geom_weights[2][4, iny]     * pwΔmpz
            Δz_ext_inv_r_pwΔmpy::T = Δz_ext_inv_r * pwΔmpy
            Δz_ext_inv_l_pwΔmpy::T = Δz_ext_inv_l * pwΔmpy

            pwwyr_pwwzr::T = pwwyr * pwwzr
            pwwyr_pwwzl::T = pwwyr * pwwzl
            pwwyl_pwwzr::T = pwwyl * pwwzr
            pwwyl_pwwzl::T = pwwyl * pwwzl

            load_weights_for_innerloop!(line_weights, pssrb, iy, iny, iz, inz,
                update_even_points, Val(idx3_is_even ? true : false), 
                pwwzr, pwwzl, pwwyr, pwwyl,
                pwwyr_pwwzr, pwwyr_pwwzl, pwwyl_pwwzr, pwwyl_pwwzl,
                pwΔmpy_pwΔmpz,
                Δy_ext_inv_r_pwΔmpz, Δy_ext_inv_l_pwΔmpz, 
                Δz_ext_inv_r_pwΔmpy, Δz_ext_inv_l_pwΔmpy,
            )

            innerloop!(line_weights, pssrb, iy, iny, iz, inz, rb_tar_idx, rb_src_idx,
                update_even_points, Val(idx3_is_even ? true : false), 
                depletion_handling, is_weighting_potential, only2d)
        end 
        for iy in 3:2:(size(pssrb.potential, 2) - 1)
            iny::Int = iy - 1

            pwwyr::T  = pssrb.geom_weights[2][1, iny]
            pwwyl::T  = pssrb.geom_weights[2][2, iny]
            pwΔmpy::T = pssrb.geom_weights[2][3, iny] 
            pwΔmpy_pwΔmpz::T = pwΔmpy * pwΔmpz
            Δy_ext_inv_r_pwΔmpz::T  = pssrb.geom_weights[2][4, iny + 1] * pwΔmpz
            Δy_ext_inv_l_pwΔmpz::T  = pssrb.geom_weights[2][4, iny]     * pwΔmpz
            Δz_ext_inv_r_pwΔmpy::T = Δz_ext_inv_r * pwΔmpy
            Δz_ext_inv_l_pwΔmpy::T = Δz_ext_inv_l * pwΔmpy

            pwwyr_pwwzr::T = pwwyr * pwwzr
            pwwyr_pwwzl::T = pwwyr * pwwzl
            pwwyl_pwwzr::T = pwwyl * pwwzr
            pwwyl_pwwzl::T = pwwyl * pwwzl

            load_weights_for_innerloop!(line_weights, pssrb, iy, iny, iz, inz,
                update_even_points, Val(idx3_is_even ? false : true), 
                pwwzr, pwwzl, pwwyr, pwwyl,
                pwwyr_pwwzr, pwwyr_pwwzl, pwwyl_pwwzr, pwwyl_pwwzl,
                pwΔmpy_pwΔmpz,
                Δy_ext_inv_r_pwΔmpz, Δy_ext_inv_l_pwΔmpz, 
                Δz_ext_inv_r_pwΔmpy, Δz_ext_inv_l_pwΔmpy,
            )

            innerloop!(line_weights, pssrb, iy, iny, iz, inz, rb_tar_idx, rb_src_idx,
                update_even_points, Val(idx3_is_even ? false : true), 
                depletion_handling, is_weighting_potential, only2d)
        end 
    end 
end

function load_weights_for_innerloop!(line_weights, pssrb::PotentialSimulationSetupRB{T, Cartesian, 3, Array{T, 3}},
        iy, iny, iz, inz, 
        update_even_points, zyi_is_even_t::Val{zyi_is_even}, 
        pwwzr, pwwzl, pwwyr, pwwyl,
        pwwyr_pwwzr, pwwyr_pwwzl, pwwyl_pwwzr, pwwyl_pwwzl,
        pwΔmpy_pwΔmpz,
        Δy_ext_inv_r_pwΔmpz, Δy_ext_inv_l_pwΔmpz, 
        Δz_ext_inv_r_pwΔmpy, Δz_ext_inv_l_pwΔmpy
        ) where {T, zyi_is_even}
    @fastmath @inbounds @simd ivdep for ix in 2:(size(pssrb.potential, 1) - 1)
        inx::Int = nidx(ix, update_even_points, zyi_is_even_t)::Int

        pwwxr::T        = pssrb.geom_weights[1][1, inx]
        pwwxl::T        = pssrb.geom_weights[1][2, inx]
        pwΔmpx::T       = pssrb.geom_weights[1][3, inx]
        Δx_ext_inv_r::T = pssrb.geom_weights[1][4, inx + 1]
        Δx_ext_inv_l::T = pssrb.geom_weights[1][4, inx]
        
        pwwxr_pwwzr::T = pwwxr * pwwzr
        pwwxl_pwwzr::T = pwwxl * pwwzr
        pwwxr_pwwzl::T = pwwxr * pwwzl
        pwwxl_pwwzl::T = pwwxl * pwwzl
        pwwxl_pwwyr::T = pwwxl * pwwyr
        pwwxr_pwwyr::T = pwwxr * pwwyr
        pwwxl_pwwyl::T = pwwxl * pwwyl
        pwwxr_pwwyl::T = pwwxr * pwwyl

        ϵ_rrr::T = pssrb.ϵ_r[ inx + 1,  iy, iz]
        ϵ_rlr::T = pssrb.ϵ_r[ inx + 1, iny, iz]
        ϵ_rrl::T = pssrb.ϵ_r[ inx + 1,  iy, inz ]
        ϵ_rll::T = pssrb.ϵ_r[ inx + 1, iny, inz ]
        ϵ_lrr::T = pssrb.ϵ_r[ inx,  iy,  iz ]
        ϵ_llr::T = pssrb.ϵ_r[ inx, iny,  iz ]
        ϵ_lrl::T = pssrb.ϵ_r[ inx,  iy, inz ] 
        ϵ_lll::T = pssrb.ϵ_r[ inx, iny, inz ] 

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

        wxr *= pwΔmpy_pwΔmpz * Δx_ext_inv_r 
        wxl *= pwΔmpy_pwΔmpz * Δx_ext_inv_l
        wyr *= Δy_ext_inv_r_pwΔmpz * pwΔmpx
        wyl *= Δy_ext_inv_l_pwΔmpz * pwΔmpx
        wzr *= Δz_ext_inv_r_pwΔmpy * pwΔmpx 
        wzl *= Δz_ext_inv_l_pwΔmpy * pwΔmpx 

        line_weights[ix-1, 1] = wxr 
        line_weights[ix-1, 2] = wxl 
        line_weights[ix-1, 3] = wyr 
        line_weights[ix-1, 4] = wyl 
        line_weights[ix-1, 5] = wzr 
        line_weights[ix-1, 6] = wzl 
    end
end



function innerloop!(line_weights, pssrb::PotentialSimulationSetupRB{T, Cartesian, 3, Array{T, 3}},
        iy, iny, iz, inz, rb_tar_idx, rb_src_idx,
        update_even_points, zyi_is_even_t::Val{zyi_is_even}, 
        depletion_handling::Val{depletion_handling_enabled},
        is_weighting_potential::Val{_is_weighting_potential}, 
        only2d::Val{only_2d}) where {T, depletion_handling_enabled, _is_weighting_potential, only_2d, zyi_is_even}
    @fastmath @inbounds @simd ivdep for ix in 2:(size(pssrb.potential, 1) - 1)
    # for ix in 2:(size(pssrb.potential, 1) - 1)
        ixr::Int = get_rbidx_right_neighbour(ix, update_even_points, zyi_is_even_t)

        wxr = line_weights[ix-1, 1] 
        wxl = line_weights[ix-1, 2] 
        wyr = line_weights[ix-1, 3] 
        wyl = line_weights[ix-1, 4] 
        wzr = line_weights[ix-1, 5] 
        wzl = line_weights[ix-1, 6] 

        vxr::T = pssrb.potential[    ixr,     iy,     iz, rb_src_idx] 
        vxl::T = pssrb.potential[ixr - 1,     iy,     iz, rb_src_idx]
        vyr::T = pssrb.potential[     ix, iy + 1,     iz, rb_src_idx]
        vyl::T = pssrb.potential[     ix,    iny,     iz, rb_src_idx]
        vzr::T = pssrb.potential[     ix,     iy, iz + 1, rb_src_idx]
        vzl::T = pssrb.potential[     ix,     iy,    inz, rb_src_idx]

        new_potential::T = _is_weighting_potential ? zero(T) : (pssrb.q_eff_imp[ix, iy, iz, rb_tar_idx] + pssrb.q_eff_fix[ix, iy, iz, rb_tar_idx])
        new_potential = muladd( wxr, vxr, new_potential)
        new_potential = muladd( wxl, vxl, new_potential)
        new_potential = muladd( wyr, vyr, new_potential)
        new_potential = muladd( wyl, vyl, new_potential)
        new_potential = muladd( wzr, vzr, new_potential)
        new_potential = muladd( wzl, vzl, new_potential)

        new_potential *= pssrb.volume_weights[ix, iy, iz, rb_tar_idx]

        old_potential::T = pssrb.potential[ix, iy, iz, rb_tar_idx]

        new_potential -= old_potential
        new_potential = muladd(new_potential, pssrb.sor_const[1], old_potential)

        if depletion_handling_enabled
            vmin::T = min(vxr, vxl, vyr, vyl, vzr, vzl)
            vmax::T = max(vxr, vxl, vyr, vyl, vzr, vzl)
            
            new_point_type = pssrb.point_types[ix, iy, iz, rb_tar_idx] 
            if new_potential < vmin || new_potential > vmax
                new_potential -= pssrb.q_eff_imp[ix, iy, iz, rb_tar_idx] * pssrb.volume_weights[ix, iy, iz, rb_tar_idx] * pssrb.sor_const[1]
                if (pssrb.point_types[ix, iy, iz, rb_tar_idx] & undepleted_bit == 0) 
                    # pssrb.point_types[ix, iy, iz, rb_tar_idx] += undepleted_bit 
                    new_point_type += undepleted_bit
                end # mark this point as undepleted
            elseif pssrb.point_types[ix, iy, iz, rb_tar_idx] & undepleted_bit > 0
                new_point_type -= undepleted_bit
            end
            pssrb.point_types[ix, iy, iz, rb_tar_idx] = new_point_type
        end 

        pssrb.potential[ix, iy, iz, rb_tar_idx]::T = ifelse(pssrb.point_types[ix, iy, iz, rb_tar_idx] & update_bit > 0, new_potential, old_potential)
    end
end
