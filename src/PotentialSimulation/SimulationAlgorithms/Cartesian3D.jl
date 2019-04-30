"""
    innerloops!(  iz::Int, rb_tar_idx::Int, rb_src_idx::Int, gw_x::Array{T, 2}, gw_y::Array{T, 2}, gw_z::Array{T, 2}, fssrb::PotentialSimulationSetupRB{T, 3, 4, :cartesian},
                                update_even_points::Val{even_points},
                                depletion_handling::Val{depletion_handling_enabled},
                                bulk_is_ptype::Val{_bulk_is_ptype}  )::Nothing where {T, even_points, depletion_handling_enabled, _bulk_is_ptype}

(Vectorized) inner loop for Cartesian coordinates. This function does all the work in the field calculation.                            
"""
@fastmath function innerloops!( iz::Int, rb_tar_idx::Int, rb_src_idx::Int, gw_x::Array{T, 2}, gw_y::Array{T, 2}, gw_z::Array{T, 2}, fssrb::PotentialSimulationSetupRB{T, 3, 4, :cartesian},
                                update_even_points::Val{even_points},
                                depletion_handling::Val{depletion_handling_enabled},
                                bulk_is_ptype::Val{_bulk_is_ptype}, 
                                is_weighting_potential::Val{_is_weighting_potential}  )::Nothing where {T, even_points, depletion_handling_enabled, _bulk_is_ptype, _is_weighting_potential}
    @inbounds begin 
        inz::Int = iz - 1 
                
        pwwzr::T        = gw_z[1, inz]
        pwwzl::T        = gw_z[2, inz]
        pwΔmpz::T       = gw_z[3, inz]
        Δz_ext_inv_r::T = gw_z[4, inz + 1]
        Δz_ext_inv_l::T = gw_z[4, inz]

        for iy in 2:(size(fssrb.potential, 2) - 1)
            iny::Int = iy - 1
            zyi_is_even::Bool = iseven(iz + iy)
            zyi_is_even_t = zyi_is_even ? Val{true}() : Val{false}()

            pwwyr::T  = gw_y[1, iny]
            pwwyl::T  = gw_y[2, iny]
            pwΔmpy::T = gw_y[3, iny] 
            pwΔmpy_pwΔmpz::T = pwΔmpy * pwΔmpz
            Δy_ext_inv_r_pwΔmpz::T  = gw_y[4, iny + 1] * pwΔmpz
            Δy_ext_inv_l_pwΔmpz::T  = gw_y[4, iny]     * pwΔmpz
            Δz_ext_inv_r_pwΔmpy::T = Δz_ext_inv_r * pwΔmpy
            Δz_ext_inv_l_pwΔmpy::T = Δz_ext_inv_l * pwΔmpy

            pwwyr_pwwzr::T = pwwyr * pwwzr
            pwwyr_pwwzl::T = pwwyr * pwwzl
            pwwyl_pwwzr::T = pwwyl * pwwzr
            pwwyl_pwwzl::T = pwwyl * pwwzl

            # for ix in 2:(size(fssrb.potential, 1) - 1)
            @fastmath @inbounds @simd ivdep for ix in 2:(size(fssrb.potential, 1) - 1)
                inx::Int = nidx(ix, update_even_points, zyi_is_even_t)::Int
                ixr::Int = ifelse( zyi_is_even, ix, even_points ? ix - 1 : ix + 1)
                ixr += ifelse(even_points, 1, 0)
                pwwxr::T        = gw_x[1, inx]
                pwwxl::T        = gw_x[2, inx]
                pwΔmpx::T       = gw_x[3, inx]
                Δx_ext_inv_r::T = gw_x[4, inx + 1]
                Δx_ext_inv_l::T = gw_x[4, inx]
                
                pwwxr_pwwzr::T = pwwxr * pwwzr
                pwwxl_pwwzr::T = pwwxl * pwwzr
                pwwxr_pwwzl::T = pwwxr * pwwzl
                pwwxl_pwwzl::T = pwwxl * pwwzl
                pwwxl_pwwyr::T = pwwxl * pwwyr
                pwwxr_pwwyr::T = pwwxr * pwwyr
                pwwxl_pwwyl::T = pwwxl * pwwyl
                pwwxr_pwwyl::T = pwwxr * pwwyl

                ϵ_rrr::T = fssrb.ϵ[ inx + 1,  iy, iz]
                ϵ_rlr::T = fssrb.ϵ[ inx + 1, iny, iz]
                ϵ_rrl::T = fssrb.ϵ[ inx + 1,  iy, inz ]
                ϵ_rll::T = fssrb.ϵ[ inx + 1, iny, inz ]
                ϵ_lrr::T = fssrb.ϵ[ inx,  iy,  iz ]
                ϵ_llr::T = fssrb.ϵ[ inx, iny,  iz ]
                ϵ_lrl::T = fssrb.ϵ[ inx,  iy, inz ] 
                ϵ_lll::T = fssrb.ϵ[ inx, iny, inz ] 

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

                vxr::T = fssrb.potential[    ixr,     iy,     iz, rb_src_idx] 
                vxl::T = fssrb.potential[ixr - 1,     iy,     iz, rb_src_idx]
                vyr::T = fssrb.potential[     ix, iy + 1,     iz, rb_src_idx]
                vyl::T = fssrb.potential[     ix,    iny,     iz, rb_src_idx]
                vzr::T = fssrb.potential[     ix,     iy, iz + 1, rb_src_idx]
                vzl::T = fssrb.potential[     ix,     iy,    inz, rb_src_idx]

                new_potential::T = _is_weighting_potential ? 0 : fssrb.ρ[ix, iy, iz, rb_tar_idx]
                new_potential = muladd( wxr, vxr, new_potential)
                new_potential = muladd( wxl, vxl, new_potential)
                new_potential = muladd( wyr, vyr, new_potential)
                new_potential = muladd( wyl, vyl, new_potential)
                new_potential = muladd( wzr, vzr, new_potential)
                new_potential = muladd( wzl, vzl, new_potential)

                new_potential *= fssrb.volume_weights[ix, iy, iz, rb_tar_idx]

                old_potential::T = fssrb.potential[ix, iy, iz, rb_tar_idx]

                new_potential -= old_potential
                new_potential = muladd(new_potential, fssrb.sor_const[1], old_potential)

                if depletion_handling_enabled
                    if _bulk_is_ptype # p-type detectors
                        if new_potential < fssrb.minimum_applied_potential
                            # new_potential = fssrb.minimum_applied_potential
                            new_potential -= fssrb.ρ[ix, iy, iz, rb_tar_idx] * fssrb.volume_weights[ix, iy, iz, rb_tar_idx] * fssrb.sor_const[1]
                            if (fssrb.pointtypes[ix, iy, iz, rb_tar_idx] & undepleted_bit == 0) fssrb.pointtypes[ix, iy, iz, rb_tar_idx] += undepleted_bit end # mark this point as undepleted
                        else
                            vmin::T = ifelse( vxr <  vxl, vxr,  vxl)
                            vmin    = ifelse( vyr < vmin, vyr, vmin)
                            vmin    = ifelse( vyl < vmin, vyl, vmin)
                            vmin    = ifelse( vzr < vmin, vzr, vmin)
                            vmin    = ifelse( vzl < vmin, vzl, vmin)
                            if new_potential <= vmin # bubble point
                                new_potential -= fssrb.ρ[ix, iy, iz, rb_tar_idx] * fssrb.volume_weights[ix, iy, iz, rb_tar_idx] * fssrb.sor_const[1]
                                if (fssrb.pointtypes[ix, iy, iz, rb_tar_idx] & undepleted_bit == 0) fssrb.pointtypes[ix, iy, iz, rb_tar_idx] += undepleted_bit end # mark this point as undepleted
                            else # normal point
                                if (fssrb.pointtypes[ix, iy, iz, rb_tar_idx] & undepleted_bit > 0) fssrb.pointtypes[ix, iy, iz, rb_tar_idx] -= undepleted_bit end # unmark this point
                            end
                        end
                    else # n-type detectors
                        if new_potential > fssrb.maximum_applied_potential
                            new_potential -= fssrb.ρ[ix, iy, iz, rb_tar_idx] * fssrb.volume_weights[ix, iy, iz, rb_tar_idx] * fssrb.sor_const[1]
                            if (fssrb.pointtypes[ix, iy, iz, rb_tar_idx] & undepleted_bit == 0) fssrb.pointtypes[ix, iy, iz, rb_tar_idx] += undepleted_bit end # mark this point as undepleted
                        else
                            vmax::T = ifelse( vxr >  vxl, vxr,  vxl)
                            vmax    = ifelse( vyr > vmax, vyr, vmax)
                            vmax    = ifelse( vyl > vmax, vyl, vmax)
                            vmax    = ifelse( vzr > vmax, vzr, vmax)
                            vmax    = ifelse( vzl > vmax, vzl, vmax)
                            if new_potential >= vmax # bubble point
                                new_potential -= fssrb.ρ[ix, iy, iz, rb_tar_idx] * fssrb.volume_weights[ix, iy, iz, rb_tar_idx] * fssrb.sor_const[1]
                                if (fssrb.pointtypes[ix, iy, iz, rb_tar_idx] & undepleted_bit == 0) fssrb.pointtypes[ix, iy, iz, rb_tar_idx] += undepleted_bit end # mark this point as undepleted
                            else # normal point -> unmark
                                if (fssrb.pointtypes[ix, iy, iz, rb_tar_idx] & undepleted_bit > 0) fssrb.pointtypes[ix, iy, iz, rb_tar_idx] -= undepleted_bit end # unmark this point
                            end
                        end
                    end
                end 

                fssrb.potential[ix, iy, iz, rb_tar_idx]::T = ifelse(fssrb.pointtypes[ix, iy, iz, rb_tar_idx] & update_bit > 0, new_potential, old_potential)
            end # x loop
        end # y loop
    end # inbounds
end

