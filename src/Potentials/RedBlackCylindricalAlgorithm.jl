@inline function update_even_inner_loop!(grid::CylindricalRedBlackGrid, pw::PrecalculatedWeightsCylindricalRedBlack,
                                        inφ::Int, iz::Int, inz_even::Int, inz_odd::Int, vzm_idx_even::Int, vzm_idx_odd::Int, iφ::Int, 
                                        pwwφr::AbstractFloat, pwwφl::AbstractFloat, pwΔmpφ::AbstractFloat,
                                        Δφ_ext_inv_r::AbstractFloat, Δφ_ext_inv_l::AbstractFloat, 
                                        pwwzr_even::AbstractFloat, pwwzr_odd::AbstractFloat, 
                                        pwwzl_even::AbstractFloat, pwwzl_odd::AbstractFloat, 
                                        pwΔmpz_even::AbstractFloat, pwΔmpz_odd::AbstractFloat, 
                                        Δz_ext_inv_l_even::AbstractFloat, Δz_ext_inv_l_odd::AbstractFloat, 
                                        Δz_ext_inv_r_even::AbstractFloat, Δz_ext_inv_r_odd::AbstractFloat,
                                        depletion_handling::Val{depletion_handling_enabled})::Nothing where {depletion_handling_enabled}
    T = eltype(grid)

    @inbounds @fastmath @simd ivdep for ir in pw.r_range
        inr = ir - 1 

        is_even = iseven(ir + iφ)
        inz = is_even ? inz_even : inz_odd
        pwrbinds = inr, inφ, iz - 1, 2 # not an extended grid -> -1
        vzm_idx = ifelse( is_even, vzm_idx_even, vzm_idx_odd)

        pt   = pw.point_types[pwrbinds...]
        vm   = pw.volume_weights[pwrbinds...]
        charge = pw.charge[pwrbinds...]

        # neighbor potential values:
        vrr = grid.potential[ ir + 1,     iφ,          iz, 1]
        vrl = grid.potential[    inr,     iφ,          iz, 1]
        vφr = grid.potential[     ir, iφ + 1,          iz, 1]
        vφl = grid.potential[     ir,    inφ,          iz, 1]
        vzr = grid.potential[     ir,     iφ, vzm_idx + 1, 1] # vzm_idx = iz or iz - 1
        vzl = grid.potential[     ir,     iφ, vzm_idx    , 1]

        pwwrr        = pw.wr[1, inr]
        pwwrl        = pw.wr[2, inr]
        pwΔmpr       = pw.wr[3, inr]
        pwmprl       = pw.wr[4, inr]
        Δr_ext_inv_l = pw.wr[5, inr]
        r_inv        = pw.wr[6, inr]
        Δmpr_squared = pw.wr[7, inr]
        pwmprr       = pw.wr[4,  ir]
        Δr_ext_inv_r = pw.wr[5,  ir]

        pwwzr        = pw.wz[1, inz]
        pwwzl        = pw.wz[2, inz]
        pwΔmpz       = pw.wz[3, inz]
        Δz_ext_inv_r = pw.wz[4, inz + 1]
        Δz_ext_inv_l = pw.wz[4, inz]

        # right weight in r: wrr
        wrr     = pw.epsilon_r[  ir,  iφ, inz + 1] * pwwφr * pwwzr
        wrr    += pw.epsilon_r[  ir, inφ, inz + 1] * pwwφl * pwwzr   
        wrr    += pw.epsilon_r[  ir,  iφ, inz ]    * pwwφr * pwwzl    
        wrr    += pw.epsilon_r[  ir, inφ, inz ]    * pwwφl * pwwzl
        # # left weight in r: wrr
        wrl     = pw.epsilon_r[ inr,  iφ, inz + 1] * pwwφr * pwwzr
        wrl    += pw.epsilon_r[ inr, inφ, inz + 1] * pwwφl * pwwzr   
        wrl    += pw.epsilon_r[ inr,  iφ, inz ]    * pwwφr * pwwzl    
        wrl    += pw.epsilon_r[ inr, inφ, inz ]    * pwwφl * pwwzl 
        # right weight in φ: wφr
        wφr     = pw.epsilon_r[ inr,  iφ, inz + 1] * pwwrl * pwwzr 
        wφr    += pw.epsilon_r[  ir,  iφ, inz + 1] * pwwrr * pwwzr  
        wφr    += pw.epsilon_r[ inr,  iφ, inz ]    * pwwrl * pwwzl    
        wφr    += pw.epsilon_r[  ir,  iφ, inz ]    * pwwrr * pwwzl 
        # left weight in φ: wφl
        wφl     = pw.epsilon_r[ inr, inφ, inz + 1] * pwwrl * pwwzr 
        wφl    += pw.epsilon_r[  ir, inφ, inz + 1] * pwwrr * pwwzr  
        wφl    += pw.epsilon_r[ inr, inφ, inz ]    * pwwrl * pwwzl    
        wφl    += pw.epsilon_r[  ir, inφ, inz ]    * pwwrr * pwwzl 
        # right weight in z: wzr
        wzr     = pw.epsilon_r[  ir,  iφ, inz + 1] * pwwrr * pwwφr  
        wzr    += pw.epsilon_r[  ir, inφ, inz + 1] * pwwrr * pwwφl     
        wzr    += pw.epsilon_r[ inr,  iφ, inz + 1] * pwwrl * pwwφr     
        wzr    += pw.epsilon_r[ inr, inφ, inz + 1] * pwwrl * pwwφl
        # left weight in z: wzr
        wzl     = pw.epsilon_r[  ir,  iφ,    inz ] * pwwrr * pwwφr 
        wzl    += pw.epsilon_r[  ir, inφ,    inz ] * pwwrr * pwwφl    
        wzl    += pw.epsilon_r[ inr,  iφ,    inz ] * pwwrl * pwwφr    
        wzl    += pw.epsilon_r[ inr, inφ,    inz ] * pwwrl * pwwφl

        wrr    *= Δr_ext_inv_r * pwmprr * pwΔmpφ * pwΔmpz
        wrl    *= Δr_ext_inv_l * pwmprl * pwΔmpφ * pwΔmpz
        wφr    *= r_inv * Δφ_ext_inv_r * pwΔmpr * pwΔmpz
        wφl    *= r_inv * Δφ_ext_inv_l * pwΔmpr * pwΔmpz
        wzr    *= Δz_ext_inv_r * pwΔmpφ * Δmpr_squared
        wzl    *= Δz_ext_inv_l * pwΔmpφ * Δmpr_squared

        new_potential = charge

        if inr == 1 wrl = 0 end

        new_potential = muladd( wrr, vrr, new_potential)
        new_potential = muladd( wrl, vrl, new_potential)
        new_potential = muladd( wφr, vφr, new_potential)
        new_potential = muladd( wφl, vφl, new_potential)
        new_potential = muladd( wzr, vzr, new_potential)
        new_potential = muladd( wzl, vzl, new_potential)

        new_potential *= vm

        old_potential::T = grid.potential[ir, iφ, iz, 2]

        new_potential -= old_potential
        new_potential = muladd(new_potential, pw.sor_const[inr], old_potential)

        if depletion_handling_enabled
            if new_potential < pw.minimum_voltage
                # new_potential = pw.minimum_voltage
                new_potential = pw.minimum_voltage
                if (pt & undepleted_bit == 0) pw.point_types[pwrbinds...] = pw.point_types[pwrbinds...] + undepleted_bit end # mark this point as undepleted
            else
                vmin = ifelse( vrr <  vrl, vrr,  vrl)
                vmin = ifelse( vφr < vmin, vφr, vmin)
                vmin = ifelse( vφl < vmin, vφl, vmin)
                vmin = ifelse( vzr < vmin, vzr, vmin)
                vmin = ifelse( vzl < vmin, vzl, vmin)
                if new_potential <= vmin # bubble point
                    new_potential -= charge * vm * pw.sor_const[inr]
                    if (pt & undepleted_bit == 0) pw.point_types[pwrbinds...] = pw.point_types[pwrbinds...] + undepleted_bit end # mark this point as undepleted
                else # normal point
                    if (pt & undepleted_bit > 0) pw.point_types[pwrbinds...] = pw.point_types[pwrbinds...] - undepleted_bit end # unmark this point
                end
            end
        end
    
        grid.potential[ir, iφ, iz, 2] = ifelse(pt & update_bit > 0, new_potential, old_potential)
    end
    return nothing
end


@inline function update_odd_inner_loop!(grid::CylindricalRedBlackGrid, pw::PrecalculatedWeightsCylindricalRedBlack,
                                        inφ::Int, iz::Int, inz_even::Int, inz_odd::Int, vzm_idx_even::Int, vzm_idx_odd::Int, iφ::Int, 
                                        pwwφr::AbstractFloat, pwwφl::AbstractFloat, pwΔmpφ::AbstractFloat,
                                        Δφ_ext_inv_r::AbstractFloat, Δφ_ext_inv_l::AbstractFloat, 
                                        pwwzr_even::AbstractFloat, pwwzr_odd::AbstractFloat, 
                                        pwwzl_even::AbstractFloat, pwwzl_odd::AbstractFloat, 
                                        pwΔmpz_even::AbstractFloat, pwΔmpz_odd::AbstractFloat, 
                                        Δz_ext_inv_l_even::AbstractFloat, Δz_ext_inv_l_odd::AbstractFloat, 
                                        Δz_ext_inv_r_even::AbstractFloat, Δz_ext_inv_r_odd::AbstractFloat,
                                        depletion_handling::Val{depletion_handling_enabled})::Nothing where {depletion_handling_enabled}
    T = eltype(grid)
    @inbounds @fastmath @simd ivdep for ir in pw.r_range
        inr = ir - 1

        is_even = iseven(ir + iφ)
        inz = is_even ? inz_even : inz_odd
        pwrbinds = inr, inφ, iz - 1, 1 # not an extended grid -> -1
        vzm_idx = ifelse( is_even, vzm_idx_even, vzm_idx_odd)

        pt = pw.point_types[pwrbinds...]
        vm = pw.volume_weights[pwrbinds...]
        charge = pw.charge[pwrbinds...]

        # neighbor potential values:
        vrr = grid.potential[ ir + 1, iφ, iz, 2 ]
        vrl = grid.potential[ ir - 1, iφ, iz, 2 ]
        vφr = grid.potential[ ir, iφ + 1, iz, 2 ]
        vφl = grid.potential[ ir, iφ - 1, iz, 2 ]
        vzr = grid.potential[ ir, iφ, vzm_idx , 2]
        vzl = grid.potential[ ir, iφ, vzm_idx - 1, 2]

        pwwrr        = pw.wr[1, inr]
        pwwrl        = pw.wr[2, inr]
        pwΔmpr       = pw.wr[3, inr]
        pwmprl       = pw.wr[4, inr]
        Δr_ext_inv_l = pw.wr[5, inr]
        r_inv        = pw.wr[6, inr]
        Δmpr_squared = pw.wr[7, inr]
        pwmprr       = pw.wr[4,  ir]
        Δr_ext_inv_r = pw.wr[5,  ir]

        pwwzr        = is_even ? pwwzr_even : pwwzr_odd 
        pwwzl        = is_even ? pwwzl_even : pwwzl_odd 
        pwΔmpz       = is_even ? pwΔmpz_even : pwΔmpz_odd 
        Δz_ext_inv_l = is_even ? Δz_ext_inv_l_even : Δz_ext_inv_l_odd 
        Δz_ext_inv_r = is_even ? Δz_ext_inv_r_even : Δz_ext_inv_r_odd 


        # right weight in r: wrr
        wrr     = pw.epsilon_r[  ir,  iφ, inz + 1] * pwwφr * pwwzr
        wrr    += pw.epsilon_r[  ir, inφ, inz + 1] * pwwφl * pwwzr   
        wrr    += pw.epsilon_r[  ir,  iφ, inz ]    * pwwφr * pwwzl    
        wrr    += pw.epsilon_r[  ir, inφ, inz ]    * pwwφl * pwwzl
        # # left weight in r: wrr
        wrl     = pw.epsilon_r[ inr,  iφ, inz + 1] * pwwφr * pwwzr
        wrl    += pw.epsilon_r[ inr, inφ, inz + 1] * pwwφl * pwwzr   
        wrl    += pw.epsilon_r[ inr,  iφ, inz ]    * pwwφr * pwwzl    
        wrl    += pw.epsilon_r[ inr, inφ, inz ]    * pwwφl * pwwzl 
        # right weight in φ: wφr
        wφr     = pw.epsilon_r[ inr,  iφ, inz + 1] * pwwrl * pwwzr 
        wφr    += pw.epsilon_r[  ir,  iφ, inz + 1] * pwwrr * pwwzr  
        wφr    += pw.epsilon_r[ inr,  iφ, inz ]    * pwwrl * pwwzl    
        wφr    += pw.epsilon_r[  ir,  iφ, inz ]    * pwwrr * pwwzl 
        # left weight in φ: wφl
        wφl     = pw.epsilon_r[ inr, inφ, inz + 1] * pwwrl * pwwzr 
        wφl    += pw.epsilon_r[  ir, inφ, inz + 1] * pwwrr * pwwzr  
        wφl    += pw.epsilon_r[ inr, inφ, inz ]    * pwwrl * pwwzl    
        wφl    += pw.epsilon_r[  ir, inφ, inz ]    * pwwrr * pwwzl 
        # right weight in z: wzr
        wzr     = pw.epsilon_r[  ir,  iφ, inz + 1] * pwwrr * pwwφr  
        wzr    += pw.epsilon_r[  ir, inφ, inz + 1] * pwwrr * pwwφl     
        wzr    += pw.epsilon_r[ inr,  iφ, inz + 1] * pwwrl * pwwφr     
        wzr    += pw.epsilon_r[ inr, inφ, inz + 1] * pwwrl * pwwφl
        # left weight in z: wzr
        wzl     = pw.epsilon_r[  ir,  iφ,    inz ] * pwwrr * pwwφr 
        wzl    += pw.epsilon_r[  ir, inφ,    inz ] * pwwrr * pwwφl    
        wzl    += pw.epsilon_r[ inr,  iφ,    inz ] * pwwrl * pwwφr    
        wzl    += pw.epsilon_r[ inr, inφ,    inz ] * pwwrl * pwwφl

        wrr    *= Δr_ext_inv_r * pwmprr * pwΔmpφ * pwΔmpz
        wrl    *= Δr_ext_inv_l * pwmprl * pwΔmpφ * pwΔmpz
        wφr    *= r_inv * Δφ_ext_inv_r * pwΔmpr * pwΔmpz
        wφl    *= r_inv * Δφ_ext_inv_l * pwΔmpr * pwΔmpz
        wzr    *= Δz_ext_inv_r * pwΔmpφ * Δmpr_squared
        wzl    *= Δz_ext_inv_l * pwΔmpφ * Δmpr_squared

        new_potential = charge

        if inr == 1 wrl = 0 end

        new_potential = muladd( wrr, vrr, new_potential)
        new_potential = muladd( wrl, vrl, new_potential)
        new_potential = muladd( wφr, vφr, new_potential)
        new_potential = muladd( wφl, vφl, new_potential)
        new_potential = muladd( wzr, vzr, new_potential)
        new_potential = muladd( wzl, vzl, new_potential)

        new_potential *= vm 

        
        old_potential = grid.potential[ir, iφ, iz, 1]               

        new_potential -= old_potential      
        new_potential = muladd(new_potential, pw.sor_const[inr], old_potential)

        if depletion_handling_enabled
            if new_potential < pw.minimum_voltage
                new_potential = pw.minimum_voltage
                if (pt & undepleted_bit == 0) pw.point_types[pwrbinds...] = pw.point_types[pwrbinds...] + undepleted_bit end # mark this point as undepleted
            else
                vmin = ifelse( vrr <  vrl, vrr,  vrl)
                vmin = ifelse( vφr < vmin, vφr, vmin)
                vmin = ifelse( vφl < vmin, vφl, vmin)
                vmin = ifelse( vzr < vmin, vzr, vmin)
                vmin = ifelse( vzl < vmin, vzl, vmin)
                if new_potential <= vmin # bubble point
                    new_potential -= charge * vm * pw.sor_const[inr]
                    if (pt & undepleted_bit == 0) pw.point_types[pwrbinds...] = pw.point_types[pwrbinds...] + undepleted_bit end # mark this point as undepleted
                else # normal point
                    if (pt & undepleted_bit > 0) pw.point_types[pwrbinds...] = pw.point_types[pwrbinds...] - undepleted_bit end # unmark this point 
                end
            end
        end

        grid.potential[ir, iφ, iz, 1] = ifelse(pt & update_bit > 0, new_potential, old_potential)

    end
    return nothing
end

function update_even_pot!(  grid::CylindricalRedBlackGrid, pw::PrecalculatedWeightsCylindricalRedBlack,
                            depletion_handling::Val{depletion_handling_enabled}, nthreads::Int)::Nothing where {depletion_handling_enabled}
    @onthreads 1:nthreads for iz in workpart(pw.z_range, 1:nthreads, Base.Threads.threadid())
    # for iz in pw.z_range
        @inbounds @fastmath begin
        # for iz in pw.z_range
            inz_evenodd = (iz - 1) * 2 
            inz_even = inz_evenodd 
            inz_odd  = inz_evenodd - 1
            vzm_idx_even = iz
            vzm_idx_odd  = iz - 1   

            pwwzr_even  = pw.wz[1, inz_even]
            pwwzr_odd   = pw.wz[1, inz_odd ]
            pwwzl_even  = pw.wz[2, inz_even]    
            pwwzl_odd   = pw.wz[2, inz_odd ]
            pwΔmpz_even = pw.wz[3, inz_even]
            pwΔmpz_odd  = pw.wz[3, inz_odd ]
            Δz_ext_inv_l_even   = pw.wz[4, inz_even]
            Δz_ext_inv_l_odd    = pw.wz[4, inz_odd ]
            Δz_ext_inv_r_even   = pw.wz[4, inz_even + 1]
            Δz_ext_inv_r_odd    = pw.wz[4, inz_odd  + 1]

            # @onthreads whichthreads for iφ in threadpartition(pw.φ_range)
            for iφ in pw.φ_range
                inφ = iφ - 1

                pwwφr        = pw.wφ[1, inφ]
                pwwφl        = pw.wφ[2, inφ]
                pwΔmpφ       = pw.wφ[3, inφ]
                Δφ_ext_inv_r = pw.wφ[4,  iφ]
                Δφ_ext_inv_l = pw.wφ[4, inφ]

                update_even_inner_loop!(grid, pw, inφ, iz, inz_even, inz_odd, vzm_idx_even, vzm_idx_odd, iφ, 
                        pwwφr, pwwφl, pwΔmpφ, Δφ_ext_inv_r, Δφ_ext_inv_l,
                        pwwzr_even, pwwzr_odd, pwwzl_even, pwwzl_odd,
                        pwΔmpz_even, pwΔmpz_odd, Δz_ext_inv_l_even, Δz_ext_inv_l_odd, Δz_ext_inv_r_even, Δz_ext_inv_r_odd,
                        depletion_handling)     
            end
        end
    end

    nothing
end

function update_odd_pot!(   grid::CylindricalRedBlackGrid, pw::PrecalculatedWeightsCylindricalRedBlack,
                            depletion_handling::Val{depletion_handling_enabled}, nthreads::Int)::Nothing where {depletion_handling_enabled}
    @onthreads 1:nthreads for iz in workpart(pw.z_range, 1:nthreads, Base.Threads.threadid())
    # for iz in pw.z_range
        @inbounds @fastmath begin 
        # for iz in pw.z_range
            inz_evenodd = (iz - 1) * 2 
            inz_even = inz_evenodd - 1
            inz_odd  = inz_evenodd
            vzm_idx_even = iz
            vzm_idx_odd  = iz + 1

            pwwzr_even  = pw.wz[1, inz_even]
            pwwzr_odd   = pw.wz[1, inz_odd ]
            pwwzl_even  = pw.wz[2, inz_even]    
            pwwzl_odd   = pw.wz[2, inz_odd ]
            pwΔmpz_even = pw.wz[3, inz_even]
            pwΔmpz_odd  = pw.wz[3, inz_odd ]
            Δz_ext_inv_l_even   = pw.wz[4, inz_even]
            Δz_ext_inv_l_odd    = pw.wz[4, inz_odd ]
            Δz_ext_inv_r_even   = pw.wz[4, inz_even + 1]
            Δz_ext_inv_r_odd    = pw.wz[4, inz_odd  + 1]

            for iφ in pw.φ_range
            # @onthreads whichthreads for iφ in threadpartition(pw.φ_range)
                inφ = iφ - 1
                
                pwwφr        = pw.wφ[1, inφ]
                pwwφl        = pw.wφ[2, inφ]
                pwΔmpφ       = pw.wφ[3, inφ]
                Δφ_ext_inv_r = pw.wφ[4,  iφ]
                Δφ_ext_inv_l = pw.wφ[4, inφ]

                update_odd_inner_loop!(grid, pw, inφ, iz, inz_even, inz_odd, vzm_idx_even, vzm_idx_odd, iφ, 
                        pwwφr, pwwφl, pwΔmpφ, Δφ_ext_inv_r, Δφ_ext_inv_l,
                        pwwzr_even, pwwzr_odd, pwwzl_even, pwwzl_odd,
                        pwΔmpz_even, pwΔmpz_odd, Δz_ext_inv_l_even, Δz_ext_inv_l_odd, Δz_ext_inv_r_even, Δz_ext_inv_r_odd,
                        depletion_handling)     
            end
        end
    end
    nothing
end

function update_pot!(   grid::CylindricalRedBlackGrid, pw::PrecalculatedWeightsCylindricalRedBlack,
                        depletion_handling::Val{depletion_handling_enabled}, nthreads::Int)::Nothing where {depletion_handling_enabled}

    update_even_pot!(grid, pw, depletion_handling, nthreads)
    
    @inbounds grid.potential[1,   :,   :, 2] = grid.potential[      3,       :,       :, 2] # reflecting boundary
    @inbounds grid.potential[end, :,   :, 2] = grid.potential[end - 2,       :,       :, 2] # reflecting boundary
    @inbounds grid.potential[:,   1,   :, 2] = grid.potential[      :, end - 1,       :, 2] # cycling boundary
    @inbounds grid.potential[:, end,   :, 2] = grid.potential[      :,       2,       :, 2] # cycling boundary
    @inbounds grid.potential[:,   :,   1, 2] = grid.potential[      :,       :,       2, 2] # reflecting boundary
    @inbounds grid.potential[:,   :, end, 2] = grid.potential[      :,       :, end - 1, 2] # reflecting boundary

    T = eltype(grid.potential)

    nφ::Int = size(pw.epsilon_r, 2) - 1
    Δmpφ::Array{T, 1} = pw.wφ[3, 1:nφ]

    @inbounds @fastmath for inz in 1:(size(pw.epsilon_r, 3) - 1)
        m::T = 0
        l::T = 0
        for inφ in 1:nφ
            inr = 1
            if iseven(inz + inφ + inr)
                rbinds = get_rb_inds(inr, inφ, inz)
                l += Δmpφ[inφ]
                m += grid.potential[rbinds..., 2] * Δmpφ[inφ]
            end
        end
        m *= inv(l)
        for inφ in 1:nφ
            inr = 1
            if iseven(inz + inφ + inr)
                rbinds = get_rb_inds(inr, inφ, inz)
                grid.potential[rbinds..., 2] = m
            end
        end
    end

    update_odd_pot!(grid, pw, depletion_handling, nthreads)

    @inbounds grid.potential[1,   :,   :, 1] = grid.potential[      3,       :,       :, 1] # reflecting boundary
    @inbounds grid.potential[end, :,   :, 1] = grid.potential[end - 2,       :,       :, 1] # reflecting boundary
    @inbounds grid.potential[:,   1,   :, 1] = grid.potential[      :, end - 1,       :, 1] # cycling boundary
    @inbounds grid.potential[:, end,   :, 1] = grid.potential[      :,       2,       :, 1] # cycling boundary
    @inbounds grid.potential[:,   :,   1, 1] = grid.potential[      :,       :,       2, 1] # reflecting boundary
    @inbounds grid.potential[:,   :, end, 1] = grid.potential[      :,       :, end - 1, 1] # reflecting boundary

    @inbounds @fastmath for inz in 1:size(pw.epsilon_r, 3)-1
        m::T = 0
        l::T = 0
        for inφ in 1:size(pw.epsilon_r, 2)-1
            inr = 1
            if isodd(inz + inφ + inr)
                rbinds = get_rb_inds(inr, inφ, inz)
                l += Δmpφ[inφ]
                m += grid.potential[rbinds..., 1] * Δmpφ[inφ]
            end
        end
        m *= inv(l)
        for inφ in 1:size(pw.epsilon_r, 2)-1
            inr = 1
            if isodd(inz + inφ + inr)
                rbinds = get_rb_inds(inr, inφ, inz)
                grid.potential[rbinds..., 1] = m
            end
        end
    end

    nothing
end

