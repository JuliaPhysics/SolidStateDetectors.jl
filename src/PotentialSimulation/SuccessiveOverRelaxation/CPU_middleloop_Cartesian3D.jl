@fastmath function middleloop!( iz::Int, rb_tar_idx::Int, rb_src_idx::Int, 
                                pssrb::PotentialSimulationSetupRB{T, Cartesian, 3, Array{T, 3}},
                                update_even_points::Val{even_points},
                                depletion_handling::Val{depletion_handling_enabled},
                                is_weighting_potential::Val{_is_weighting_potential},
                                only2d::Val{only_2d},
                                idx3iseven::Val{idx3_is_even})::Nothing where {T, even_points, depletion_handling_enabled, _is_weighting_potential, only_2d, idx3_is_even}
    @inbounds begin 
        inz::Int = iz - 1 
                
        pwwzr::T        = pssrb.geom_weights[1][1, inz]
        pwwzl::T        = pssrb.geom_weights[1][2, inz]
        pwΔmpz::T       = pssrb.geom_weights[1][3, inz]
        Δz_ext_inv_r::T = pssrb.geom_weights[1][4, inz + 1]
        Δz_ext_inv_l::T = pssrb.geom_weights[1][4, inz]

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
                Δz_ext_inv_r_pwΔmpy, Δz_ext_inv_l_pwΔmpy,
                Δy_ext_inv_r_pwΔmpz, Δy_ext_inv_l_pwΔmpz, 
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
                Δz_ext_inv_r_pwΔmpy, Δz_ext_inv_l_pwΔmpy,
                Δy_ext_inv_r_pwΔmpz, Δy_ext_inv_l_pwΔmpz, 
            )

            innerloop!(line_weights, pssrb, iy, iny, iz, inz, rb_tar_idx, rb_src_idx,
                update_even_points, Val(idx3_is_even ? false : true), 
                depletion_handling, is_weighting_potential, only2d)
        end 
    end 
end