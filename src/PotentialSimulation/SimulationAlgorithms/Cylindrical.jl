"""
    update!(fssrb::PotentialSimulationSetupRB{T, 3, 4, S}, RBT::DataType)::Nothing

Loop over `even` grid points. A point is `even` if the sum of its cartesian indicies (of the not extended grid) is even.
Even points get the red black index (rbi) = 2. ( -> rbpotential[ inds..., rbi ]).
"""
@fastmath function update!( fssrb::PotentialSimulationSetupRB{T, 3, 4, S}, use_nthreads::Int,
                            update_even_points::Val{even_points},
                            depletion_handling::Val{depletion_handling_enabled},
                            bulk_is_ptype::Val{_bulk_is_ptype},
                            is_weighting_potential::Val{_is_weighting_potential})::Nothing where {T, S, even_points, depletion_handling_enabled, _bulk_is_ptype, _is_weighting_potential}
    @inbounds begin 
        rb_tar_idx::Int, rb_src_idx::Int = even_points ? (rb_even::Int, rb_odd::Int) : (rb_odd::Int,rb_even::Int) 

        gw1::Array{T, 2} = fssrb.geom_weights[1].weights
        gw2::Array{T, 2} = fssrb.geom_weights[2].weights
        gw3::Array{T, 2} = fssrb.geom_weights[3].weights

        # for ir in 2:(size(fssrb.potential, 3) - 1)
        @onthreads 1:use_nthreads for idx1 in workpart(2:(size(fssrb.potential, 3) - 1), 1:use_nthreads, Base.Threads.threadid())
            innerloops!( idx1, rb_tar_idx, rb_src_idx, gw1, gw2, gw3, fssrb, update_even_points, depletion_handling, bulk_is_ptype, is_weighting_potential)
        end 
    end 
    nothing
end

"""
    innerloops!(  ir::Int, rb_tar_idx::Int, rb_src_idx::Int, gw_r::Array{T, 2}, gw_θ::Array{T, 2}, gw_z::Array{T, 2}, fssrb::PotentialSimulationSetupRB{T, 3, 4, :Cylindrical},
                                update_even_points::Val{even_points},
                                depletion_handling::Val{depletion_handling_enabled},
                                bulk_is_ptype::Val{_bulk_is_ptype}  )::Nothing where {T, even_points, depletion_handling_enabled, _bulk_is_ptype}

(Vectorized) inner loop for Cylindrical coordinates. This function does all the work in the fied calculation.                            
"""
@fastmath function innerloops!( ir::Int, rb_tar_idx::Int, rb_src_idx::Int, gw_r::Array{T, 2}, gw_θ::Array{T, 2}, gw_z::Array{T, 2}, fssrb::PotentialSimulationSetupRB{T, 3, 4, :Cylindrical},
                                update_even_points::Val{even_points},
                                depletion_handling::Val{depletion_handling_enabled},
                                bulk_is_ptype::Val{_bulk_is_ptype}, 
                                is_weighting_potential::Val{_is_weighting_potential}  )::Nothing where {T, even_points, depletion_handling_enabled, _bulk_is_ptype, _is_weighting_potential}
    @inbounds begin 
        inr::Int = ir - 1 
                
        pwwrr::T               = gw_r[1, inr]
        pwwrl::T               = gw_r[2, inr]
        r_inv_pwΔmpr::T        = gw_r[3, inr]
        Δr_ext_inv_r_pwmprr::T = gw_r[4, inr] 
        Δr_ext_inv_l_pwmprl::T = gw_r[5, inr] 
        Δmpr_squared::T        = gw_r[6, inr]  

        for iθ in 2:(size(fssrb.potential, 2) - 1)
            # θ loop
            inθ::Int = iθ - 1
            pwwθr::T        = gw_θ[1, inθ]
            pwwθl::T        = gw_θ[2, inθ]
            pwΔmpθ::T       = gw_θ[3, inθ]
            Δθ_ext_inv_r::T = gw_θ[4,  iθ]
            Δθ_ext_inv_l::T = gw_θ[4, inθ]

            if inr == 1
                pwwθr = 0.5f0
                pwwθl = 0.5f0
                pwΔmpθ = 2π
                Δθ_ext_inv_r = 0.15915494f0
                Δθ_ext_inv_l = 0.15915494f0
            end
            
            rθi_is_even::Bool = iseven(ir + iθ)
            rθi_is_even_t = iseven(ir + iθ) ? Val{true}() : Val{false}()

            pwwrr_pwwθr::T = pwwrr * pwwθr
            pwwrr_pwwθl::T = pwwrr * pwwθl
            pwwrl_pwwθr::T = pwwrl * pwwθr
            pwwrl_pwwθl::T = pwwrl * pwwθl

            Δr_ext_inv_r_pwmprr_pwΔmpθ::T = Δr_ext_inv_r_pwmprr * pwΔmpθ
            Δr_ext_inv_l_pwmprl_pwΔmpθ::T = Δr_ext_inv_l_pwmprl * pwΔmpθ
            pwΔmpθ_Δmpr_squared::T = pwΔmpθ * Δmpr_squared
            r_inv_pwΔmpr_Δθ_ext_inv_r::T = r_inv_pwΔmpr * Δθ_ext_inv_r
            r_inv_pwΔmpr_Δθ_ext_inv_l::T = r_inv_pwΔmpr * Δθ_ext_inv_l

            @fastmath @inbounds @simd ivdep for iz in 2:(size(fssrb.potential, 1) - 1)
                inz::Int = nidx(iz, update_even_points, rθi_is_even_t)::Int
                # izr::Int = get_rbidx_right_neighbour(iz, update_even_points, rθi_is_even)::Int # this is somehow slower than the two lines below
                izr::Int = ifelse( rθi_is_even, iz, even_points ? iz - 1 : iz + 1)
                izr += ifelse(even_points, 1, 0)

                pwwzr::T        = gw_z[1, inz]
                pwwzl::T        = gw_z[2, inz]
                pwΔmpz::T       = gw_z[3, inz]
                Δz_ext_inv_r::T = gw_z[4, inz + 1]
                Δz_ext_inv_l::T = gw_z[4, inz]

                ϵ_rrr::T = fssrb.ϵ[  ir,  iθ, inz + 1]
                ϵ_rlr::T = fssrb.ϵ[  ir, inθ, inz + 1]
                ϵ_rrl::T = fssrb.ϵ[  ir,  iθ, inz ]
                ϵ_rll::T = fssrb.ϵ[  ir, inθ, inz ]
                ϵ_lrr::T = fssrb.ϵ[ inr,  iθ, inz + 1]
                ϵ_llr::T = fssrb.ϵ[ inr, inθ, inz + 1]
                ϵ_lrl::T = fssrb.ϵ[ inr,  iθ, inz ] 
                ϵ_lll::T = fssrb.ϵ[ inr, inθ, inz ] 

                vrr::T = fssrb.potential[     iz,     iθ, ir + 1, rb_src_idx]
                vrl::T = fssrb.potential[     iz,     iθ,    inr, rb_src_idx]
                vθr::T = fssrb.potential[     iz, iθ + 1,     ir, rb_src_idx]
                vθl::T = fssrb.potential[     iz,    inθ,     ir, rb_src_idx]
                vzr::T = fssrb.potential[    izr,     iθ,     ir, rb_src_idx] 
                vzl::T = fssrb.potential[izr - 1,     iθ,     ir, rb_src_idx]

                pwwθr_pwwzr::T = pwwθr * pwwzr
                pwwθl_pwwzr::T = pwwθl * pwwzr
                pwwθr_pwwzl::T = pwwθr * pwwzl
                pwwθl_pwwzl::T = pwwθl * pwwzl
                pwwrl_pwwzr::T = pwwrl * pwwzr
                pwwrr_pwwzr::T = pwwrr * pwwzr
                pwwrl_pwwzl::T = pwwrl * pwwzl
                pwwrr_pwwzl::T = pwwrr * pwwzl

                # right weight in r: wrr
                wrr::T = ϵ_rrr * pwwθr_pwwzr
                wrr    = muladd(ϵ_rlr, pwwθl_pwwzr, wrr)   
                wrr    = muladd(ϵ_rrl, pwwθr_pwwzl, wrr)    
                wrr    = muladd(ϵ_rll, pwwθl_pwwzl, wrr)
                # left weight in r: wrr
                wrl::T = ϵ_lrr * pwwθr_pwwzr
                wrl    = muladd(ϵ_llr, pwwθl_pwwzr, wrl)   
                wrl    = muladd(ϵ_lrl, pwwθr_pwwzl, wrl)    
                wrl    = muladd(ϵ_lll, pwwθl_pwwzl, wrl) 
                # right weight in θ: wθr
                wθr::T = ϵ_lrr * pwwrl_pwwzr 
                wθr    = muladd(ϵ_rrr, pwwrr_pwwzr, wθr)  
                wθr    = muladd(ϵ_lrl, pwwrl_pwwzl, wθr)    
                wθr    = muladd(ϵ_rrl, pwwrr_pwwzl, wθr) 
                # left weight in θ: wθl
                wθl::T = ϵ_llr * pwwrl_pwwzr 
                wθl    = muladd(ϵ_rlr, pwwrr_pwwzr, wθl)  
                wθl    = muladd(ϵ_lll, pwwrl_pwwzl, wθl)    
                wθl    = muladd(ϵ_rll, pwwrr_pwwzl, wθl) 
                # right weight in z: wzr
                wzr::T = ϵ_rrr * pwwrr_pwwθr  
                wzr    = muladd(ϵ_rlr, pwwrr_pwwθl, wzr)     
                wzr    = muladd(ϵ_lrr, pwwrl_pwwθr, wzr)     
                wzr    = muladd(ϵ_llr, pwwrl_pwwθl, wzr)
                # left weight in z: wzr
                wzl::T = ϵ_rrl * pwwrr_pwwθr 
                wzl    = muladd(ϵ_rll, pwwrr_pwwθl, wzl)    
                wzl    = muladd(ϵ_lrl, pwwrl_pwwθr, wzl)    
                wzl    = muladd(ϵ_lll, pwwrl_pwwθl, wzl)

                wrr *= Δr_ext_inv_r_pwmprr_pwΔmpθ * pwΔmpz
                wrl *= Δr_ext_inv_l_pwmprl_pwΔmpθ * pwΔmpz
                wθr *= r_inv_pwΔmpr_Δθ_ext_inv_r * pwΔmpz
                wθl *= r_inv_pwΔmpr_Δθ_ext_inv_l * pwΔmpz
                wzr *= Δz_ext_inv_r * pwΔmpθ_Δmpr_squared
                wzl *= Δz_ext_inv_l * pwΔmpθ_Δmpr_squared
            
                new_potential::T = _is_weighting_potential ? 0 : fssrb.ρ[iz, iθ, ir, rb_tar_idx]
                new_potential = muladd( wrr, vrr, new_potential)
                new_potential = muladd( wrl, vrl, new_potential)
                new_potential = muladd( wθr, vθr, new_potential)
                new_potential = muladd( wθl, vθl, new_potential)
                new_potential = muladd( wzr, vzr, new_potential)
                new_potential = muladd( wzl, vzl, new_potential)

                new_potential *= fssrb.volume_weights[iz, iθ, ir, rb_tar_idx]

                old_potential::T = fssrb.potential[iz, iθ, ir, rb_tar_idx]

                new_potential -= old_potential
                new_potential = muladd(new_potential, fssrb.sor_const[inr], old_potential)

                if depletion_handling_enabled
                    if inr == 1 vrl = vrr end
                    if _bulk_is_ptype # p-type detectors
                        if new_potential < fssrb.minimum_applied_potential
                            # new_potential = fssrb.minimum_applied_potential
                            new_potential -= fssrb.ρ[iz, iθ, ir, rb_tar_idx] * fssrb.volume_weights[iz, iθ, ir, rb_tar_idx] * fssrb.sor_const[inr]
                            if (fssrb.pointtypes[iz, iθ, ir, rb_tar_idx] & undepleted_bit == 0) fssrb.pointtypes[iz, iθ, ir, rb_tar_idx] += undepleted_bit end # mark this point as undepleted
                        else
                            vmin::T = ifelse( vrr <  vrl, vrr,  vrl)
                            vmin    = ifelse( vθr < vmin, vθr, vmin)
                            vmin    = ifelse( vθl < vmin, vθl, vmin)
                            vmin    = ifelse( vzr < vmin, vzr, vmin)
                            vmin    = ifelse( vzl < vmin, vzl, vmin)
                            if new_potential <= vmin # bubble point
                                new_potential -= fssrb.ρ[iz, iθ, ir, rb_tar_idx] * fssrb.volume_weights[iz, iθ, ir, rb_tar_idx] * fssrb.sor_const[inr]
                                if (fssrb.pointtypes[iz, iθ, ir, rb_tar_idx] & undepleted_bit == 0) fssrb.pointtypes[iz, iθ, ir, rb_tar_idx] += undepleted_bit end # mark this point as undepleted
                            else # normal point
                                if (fssrb.pointtypes[iz, iθ, ir, rb_tar_idx] & undepleted_bit > 0) fssrb.pointtypes[iz, iθ, ir, rb_tar_idx] -= undepleted_bit end # unmark this point
                            end
                        end
                    else # n-type detectors
                        if new_potential > fssrb.maximum_applied_potential
                            # new_potential = fssrb.maximum_applied_potential
                            new_potential -= fssrb.ρ[iz, iθ, ir, rb_tar_idx] * fssrb.volume_weights[iz, iθ, ir, rb_tar_idx] * fssrb.sor_const[inr]
                            if (fssrb.pointtypes[iz, iθ, ir, rb_tar_idx] & undepleted_bit == 0) fssrb.pointtypes[iz, iθ, ir, rb_tar_idx] += undepleted_bit end # mark this point as undepleted
                        else
                            vmax::T = ifelse( vrr >  vrl, vrr,  vrl)
                            vmax    = ifelse( vθr > vmax, vθr, vmax)
                            vmax    = ifelse( vθl > vmax, vθl, vmax)
                            vmax    = ifelse( vzr > vmax, vzr, vmax)
                            vmax    = ifelse( vzl > vmax, vzl, vmax)
                            if new_potential >= vmax # bubble point
                                new_potential -= fssrb.ρ[iz, iθ, ir, rb_tar_idx] * fssrb.volume_weights[iz, iθ, ir, rb_tar_idx] * fssrb.sor_const[inr]
                                if (fssrb.pointtypes[iz, iθ, ir, rb_tar_idx] & undepleted_bit == 0) fssrb.pointtypes[iz, iθ, ir, rb_tar_idx] += undepleted_bit end # mark this point as undepleted
                            else # normal point -> unmark
                                if (fssrb.pointtypes[iz, iθ, ir, rb_tar_idx] & undepleted_bit > 0) fssrb.pointtypes[iz, iθ, ir, rb_tar_idx] -= undepleted_bit end # unmark this point
                            end
                        end
                    end
                end

                fssrb.potential[iz, iθ, ir, rb_tar_idx]::T = ifelse(fssrb.pointtypes[iz, iθ, ir, rb_tar_idx] & update_bit > 0, new_potential, old_potential)
            end # z loop
        end # θ loop
    end # inbounds
end

function apply_boundary_conditions_on_θ_axis!(  rbpot::Array{T, 4}, iz::Int, ir::Int, rbi::Int, ax::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :open, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[iz,   1, ir, rbi] = rbpot[ iz, end - 1, ir, rbi] # cycling boundary
    rbpot[iz, end, ir, rbi] = rbpot[ iz,       2, ir, rbi] # cycling boundary
    nothing
end
function apply_boundary_conditions_on_θ_axis!(  rbpot::Array{T, 4}, iz::Int, ir::Int, rbi::Int, ax::DiscreteAxis{T, :periodic, :periodic}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[iz,   1, ir, rbi] = rbpot[ iz, end - 1, ir, rbi] # cycling boundary
    rbpot[iz, end, ir, rbi] = rbpot[ iz,       2, ir, rbi] # cycling boundary
    nothing
end
function apply_boundary_conditions_on_r_axis!(  rbpot::Array{T, 4}, iz::Int, iθ::Int, rbi::Int, ax::DiscreteAxis{T, :r0, :infinite}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[iz, iθ, end, rbi] = grid_boundary_factors[2] * rbpot[iz, iθ, end - 2, rbi] # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_r_axis!(  rbpot::Array{T, 4}, iz::Int, iθ::Int, rbi::Int, ax::DiscreteAxis{T, :r0, :reflecting}, int::Interval{:closed, :closed, T},
                                                grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[iz, iθ, end, rbi] = rbpot[iz, iθ, end - 2, rbi] # infinity boundaries in r
    nothing
end
function apply_boundary_conditions_on_cyl_z_axis!(  rbpot::Array{T, 4}, ir::Int, iθ::Int, rbi::Int, ax::DiscreteAxis{T, :infinite, :infinite}, int::Interval{:closed, :closed, T},
                                                    grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, iθ, ir, rbi] = grid_boundary_factors[1] * rbpot[ 2, iθ, ir, rbi]  # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    rbpot[end, iθ, ir, rbi] = grid_boundary_factors[2] * rbpot[ end - 1, iθ, ir, rbi]    # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    nothing
end
function apply_boundary_conditions_on_cyl_z_axis!(  rbpot::Array{T, 4}, ir::Int, iθ::Int, rbi::Int, ax::DiscreteAxis{T, :reflecting, :reflecting}, int::Interval{:closed, :closed, T},
                                                    grid_boundary_factors::NTuple{2, T})::Nothing where {T}
    rbpot[  1, iθ, ir, rbi] = rbpot[ 2, iθ, ir, rbi]  # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    rbpot[end, iθ, ir, rbi] = rbpot[ end - 1, iθ, ir, rbi]    # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    nothing
end


function apply_boundary_conditions!(fssrb::PotentialSimulationSetupRB{T, N1, N2, :Cylindrical}, update_even_points::Val{even_points}, only2d::Val{only_2d}) where {T, N1, N2, even_points, only_2d}
    rbi::Int = even_points ? rb_even::Int : rb_odd::Int
    if only_2d
        iθ::Int = 2
        @inbounds for iz in axes(fssrb.potential, 1)
            for ir in axes(fssrb.potential, 3)
                apply_boundary_conditions_on_θ_axis!( fssrb.potential, iz, ir, rbi, fssrb.grid.axes[2], fssrb.grid.axes[2].interval, fssrb.grid_boundary_factors[2])
            end
            apply_boundary_conditions_on_r_axis!( fssrb.potential, iz, iθ, rbi, fssrb.grid.axes[1], fssrb.grid.axes[1].interval, fssrb.grid_boundary_factors[1])
        end
        @inbounds for ir in axes(fssrb.potential, 3)
            apply_boundary_conditions_on_cyl_z_axis!( fssrb.potential, ir, iθ, rbi, fssrb.grid.axes[3], fssrb.grid.axes[3].interval, fssrb.grid_boundary_factors[3])
        end
    else
        @inbounds for iz in axes(fssrb.potential, 1)
            for iθ in axes(fssrb.potential, 2)
                apply_boundary_conditions_on_r_axis!( fssrb.potential, iz, iθ, rbi, fssrb.grid.axes[1], fssrb.grid.axes[1].interval, fssrb.grid_boundary_factors[1])
            end
            for ir in axes(fssrb.potential, 3)
                apply_boundary_conditions_on_θ_axis!( fssrb.potential, iz, ir, rbi, fssrb.grid.axes[2], fssrb.grid.axes[2].interval, fssrb.grid_boundary_factors[2])
            end
        end
        @inbounds for iθ in axes(fssrb.potential, 2) # z boundaries
            for ir in axes(fssrb.potential, 3)
                apply_boundary_conditions_on_cyl_z_axis!( fssrb.potential, ir, iθ, rbi, fssrb.grid.axes[3], fssrb.grid.axes[3].interval, fssrb.grid_boundary_factors[3])
            end
        end
        begin # r = 0 handling
            nθ::Int = size(fssrb.potential, 2) - 1
            gw_θ::Array{T, 2} = fssrb.geom_weights[2].weights
            @inbounds for inz in 1:(size(fssrb.ϵ, 3) - 1)
                m::T = 0
                l::T = 0
                for inθ in 1:nθ
                    if even_points ? isodd(inz + inθ) : iseven(inz + inθ)
                        l += gw_θ[3, inθ] 
                        m += fssrb.potential[rbidx(inz), inθ + 1, 2, rbi] * gw_θ[3, inθ] 
                    end
                end
                m *= inv(l)
                for inθ in 1:nθ
                    if even_points ? isodd(inz + inθ) : iseven(inz + inθ)
                        fssrb.potential[rbidx(inz), inθ + 1, 2, rbi]::T = m
                    end
                end
            end
        end
    end
    nothing
end


function update!(   fssrb::PotentialSimulationSetupRB{T}; use_nthreads::Int = Base.Threads.nthreads(), 
                    depletion_handling::Val{depletion_handling_enabled} = Val{false}(), only2d::Val{only_2d} = Val{false}(),
                    is_weighting_potential::Val{_is_weighting_potential} = Val{false}())::Nothing where {T, depletion_handling_enabled, only_2d, _is_weighting_potential}
    update!(fssrb, use_nthreads, Val{true}(), depletion_handling, Val{fssrb.bulk_is_ptype}(), is_weighting_potential)
    apply_boundary_conditions!(fssrb, Val{true}(), only2d)
    update!(fssrb, use_nthreads, Val{false}(), depletion_handling, Val{fssrb.bulk_is_ptype}(), is_weighting_potential)
    apply_boundary_conditions!(fssrb, Val{false}(), only2d)
    nothing
end
   