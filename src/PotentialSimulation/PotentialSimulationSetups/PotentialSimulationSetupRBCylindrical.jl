function PotentialSimulationSetupRB(ssd::SolidStateDetector{T}, grid::Grid{T, 3, :Cylindrical} = Grid(ssd), 
                potential_array::Union{Missing, Array{T, 3}} = missing; weighting_potential_channel_idx::Union{Missing, Int} = missing)::PotentialSimulationSetupRB{T} where {T} 
    r0_handling::Bool = typeof(grid.axes[1]).parameters[2] == :r0 
    only_2d::Bool = length(grid.axes[2]) == 1 ? true : false
    @assert grid.axes[1][1] == 0 "Something is wrong. R-axis has `:r0`-boundary handling but first tick is $(axr[1]) and not 0."

    is_weighting_potential::Bool = !ismissing(weighting_potential_channel_idx)

    @inbounds begin
        begin # Geometrical weights of the Axes 
            nr, nθ, nz = size(grid)

            # R-axis
            axr::Vector{T} = collect(grid.axes[1]) # real grid points/ticks -> the potential at these ticks are going to be calculated
            Δr::Vector{T} = diff(axr) # difference between real ticks
            r_inv::Vector{T} = inv.(axr)
            if r0_handling r_inv[1] = inv(axr[2] * 0.5) end
            r_ext::Vector{T} = get_extended_ticks(grid.axes[1])
            Δr_ext::Vector{T} = diff(r_ext)
            Δr_ext_inv::Vector{T} = inv.(Δr_ext)
            mpr::Vector{T} = midpoints(r_ext)
            mpr_inv::Vector{T} = inv.(mpr)
            Δmpr::Vector{T} = diff(mpr)
            Δmpr_inv::Vector{T} = inv.(Δmpr)
            Δmpr_squared::Vector{T} = T(0.5) .* ((mpr[2:end].^2) .- (mpr[1:end-1].^2))
            if r0_handling 
                Δmpr_squared[1] = T(0.5) * (mpr[2]^2)
                Δr_ext_inv[1] = 0 # -> left weight @r=0 becomes 0 through this
            end
            Δhmprr::Vector{T} = mpr[2:end] - axr # distances between midpoints and real grid points (half distances -> h), needed for weights wrr, wrl, ... & dV (volume element)
            Δhmprl::Vector{T} = axr - mpr[1:end - 1]
            wrr::Vector{T} = Δmpr_inv .* Δhmprr # weights for epislon_r adding
            wrl::Vector{T} = Δmpr_inv .* Δhmprl 
            if r0_handling wrl[1]::T = 1 - wrr[1] end
            wr::Array{T, 2} = zeros(T, 6, length(wrr))
            wr[1, :] = wrr
            wr[2, :] = wrl
            wr[3, :] = r_inv .* Δmpr 
            wr[4, :] = Δr_ext_inv[2:end] .* mpr[2:end]
            wr[5, :] = Δr_ext_inv[1:length(wrr)] .* mpr[1:length(wrr)]
            wr[6, :] = Δmpr_squared
            gw_r::GeometricalRadialAxisWeights{T} = GeometricalRadialAxisWeights{T}( wr ) # Weights needed for Field Simulation loop
            r_ext_mid::T = (r_ext[end - 1] - r_ext[2]) / 2
            grid_boundary_factor_r_left::T = abs((r_ext[2] - r_ext_mid) / (r_ext[1] - r_ext_mid))
            grid_boundary_factor_r_right::T = r0_handling ? abs(r_ext[end - 1] / r_ext[end]) : abs((r_ext[end - 1] - r_ext_mid) / (r_ext[end] - r_ext_mid))
            
            # θ-axis
            axθ::Vector{T} = collect(grid.axes[2]) # real grid points/ticks -> the potential at these ticks are going to be calculated
            cyclic::T = grid.axes[2].interval.right
            Δθ::Vector{T} = diff(axθ)  # difference between real ticks
            θ_ext::Vector{T} = get_extended_ticks(grid.axes[2])
            # θ_ext::Vector{T} = if cyclic > 0
            #     get_extended_ticks(grid.axes[2])
            # else
            #     T[-2π, 0, 2π]
            # end
            Δθ_ext::Vector{T} = diff(θ_ext)
            Δθ_ext_inv::Vector{T} = inv.(Δθ_ext)
            mpθ::Vector{T} = midpoints(θ_ext)
            mpθ_inv::Vector{T} = inv.(mpθ)
            Δmpθ::Vector{T} = diff(mpθ)
            Δmpθ_inv::Vector{T} = inv.(Δmpθ)
            Δhmpθr::Vector{T} = mpθ[2:end] - axθ # distances between midpoints and real grid points (half distances -> h), needed for weights wrr, wrl, ... & dV (volume element)
            Δhmpθl::Vector{T} = axθ - mpθ[1:end - 1]
            wθr::Vector{T} = Δmpθ_inv .* Δhmpθr # weights for epislon_r adding
            wθl::Vector{T} = Δmpθ_inv .* Δhmpθl
            wθ::Array{T, 2} = zeros(T, 4, length(wθr) + 1)
            wθ[1, 1:length(wθr)] = wθr
            wθ[2, 1:length(wθr)] = wθl
            wθ[3, 1:length(wθr)] = Δmpθ
            wθ[4, :] = Δθ_ext_inv
            gw_θ::GeometricalAzimutalAxisWeights{T} = GeometricalAzimutalAxisWeights{T}( wθ ) # Weights needed for Field Simulation loop
            θ_ext_mid::T = (θ_ext[end - 1] - θ_ext[2]) / 2
            grid_boundary_factor_θ_right::T = abs((θ_ext[end - 1] - θ_ext_mid) / (θ_ext[end] - θ_ext_mid))
            grid_boundary_factor_θ_left::T = abs((θ_ext[2] - θ_ext_mid) / (θ_ext[1] - θ_ext_mid))
            
            # Z-axis
            axz::Vector{T} = collect(grid.axes[3]) # real grid points/ticks -> the potential at these ticks are going to be calculated
            Δz::Vector{T} = diff(axz) # difference between real ticks
            z_ext::Vector{T} = get_extended_ticks(grid.axes[3])
            Δz_ext::Vector{T} = diff(z_ext)
            Δz_ext_inv::Vector{T} = inv.(Δz_ext)
            mpz::Vector{T} = midpoints(z_ext)
            mpz_inv::Vector{T} = inv.(mpz)
            Δmpz::Vector{T} = diff(mpz)
            Δmpz_inv::Vector{T} = inv.(Δmpz)
            Δhmpzr::Vector{T} = mpz[2:end] - axz # distances between midpoints and real grid points (half distances -> h), needed for weights wrr, wrl, ... & dV (volume element)
            Δhmpzl::Vector{T} = axz - mpz[1:end - 1]
            wzr::Vector{T} = Δmpz_inv .* Δhmpzr # weights for epislon_r adding
            wzl::Vector{T} = Δmpz_inv .* Δhmpzl
            wz::Array{T, 2} = zeros(T, 4, length(wzr) + 1)
            wz[1, 1:length(wzr)] = wzr
            wz[2, 1:length(wzr)] = wzl
            wz[3, 1:length(wzr)] = Δmpz
            wz[4, :] = Δz_ext_inv
            gw_z::GeometricalCartesianAxisWeights{T} = GeometricalCartesianAxisWeights{T}( wz ) # Weights needed for Field Simulation loop       
            z_ext_mid::T = (z_ext[end - 1] - z_ext[2]) / 2
            grid_boundary_factor_z_right::T = abs((z_ext[end - 1] - z_ext_mid) / (z_ext[end] - z_ext_mid))
            grid_boundary_factor_z_left::T = abs((z_ext[2] - z_ext_mid) / (z_ext[1] - z_ext_mid))

            geom_weights::NTuple{3, AbstractGeometricalAxisWeights{T}} = (gw_r, gw_θ, gw_z) # Weights needed for Field Simulation loop
            grid_boundary_factors::NTuple{3, NTuple{2, T}} = ((grid_boundary_factor_r_left, grid_boundary_factor_r_right), (grid_boundary_factor_θ_left, grid_boundary_factor_θ_right), (grid_boundary_factor_z_left, grid_boundary_factor_z_right))
        end
        
        detector_material_ϵ_r::T = ssd.material_detector.ϵ_r
        environment_material_ϵ_r::T = ssd.material_environment.ϵ_r
        
        bulk_is_ptype::Bool = ssd.bulk_type == :ptype ? true : false
        minimum_applied_potential::T = minimum(ssd.segment_bias_voltages) 
        maximum_applied_potential::T = maximum(ssd.segment_bias_voltages) 
        bias_voltage::T = maximum_applied_potential - minimum_applied_potential
        depletion_handling_potential_limit::T = -bias_voltage
        sor_consts::Vector{T} = [1, 1]
        sor_slope = (sor_consts[2] .- sor_consts[1]) / (nr - 1 )
        sor_const::Vector{T} = T[ sor_consts[1] + (i - 1) * sor_slope for i in 1:nr]
        
        ϵ = Array{T, 3}(undef,    length(mpr), length(mpθ), length(mpz))
        ρ_tmp = Array{T, 3}(undef, length(mpr), length(mpθ), length(mpz))
        for iz in 1:size(ϵ, 3)
            pos_z::T = mpz[iz]
            for iθ in 1:size(ϵ, 2)
                pos_θ::T = mpθ[iθ]

                ir = 1
                pos_r::T = mpr[ir]
                if axr[1] == 0
                    pos_r = axr[2] * 0.5
                end
                pt::Cylindrical{T} = Cylindrical{T}(pos_r, pos_θ, pos_z)
                ρ_tmp[ir, iθ, iz]::T, ϵ[ir, iθ, iz]::T = get_ρ_and_ϵ(pt, ssd)
                
                for ir in 2:size(ϵ, 1)
                    pos_r = mpr[ir]
                    pt = Cylindrical{T}(pos_r, pos_θ, pos_z)
                    ρ_tmp[ir, iθ, iz]::T, ϵ[ir, iθ, iz]::T = get_ρ_and_ϵ(pt, ssd)
                end
            end
        end
        ϵ0_inv::T = inv(ϵ0)
        ρ_tmp *= ϵ0_inv

        volume_weights::Array{T, 4} = RBExtBy2Array(T, grid)
        ρ::Array{T, 4} = RBExtBy2Array(T, grid)
        for iz in range(2, stop = length(z_ext) - 1)
            inz::Int = iz - 1
            irbz::Int = rbidx(inz)
            pos_z::T = z_ext[iz]
            for iθ in range(2, stop = length(θ_ext) - 1)
                inθ::Int = iθ - 1
                pos_θ::T = θ_ext[iθ]
                for ir in range(2, stop = length(r_ext) - 1)
                    inr::Int = ir - 1;
                    pos_r::T = r_ext[ir]

                    rbi::Int = iseven(inr + inθ + inz) ? rb_even::Int : rb_odd::Int
                    # rbinds = irbz, iθ, ir, rbi

                    ρ_cell::T = 0
                    if !is_weighting_potential
                        if inr > 1
                            ρ_cell += ρ_tmp[ ir,  iθ,  iz] * wzr[inz] * wrr[inr] * wθr[inθ]
                            ρ_cell += ρ_tmp[ ir,  iθ, inz] * wzl[inz] * wrr[inr] * wθr[inθ]
                            ρ_cell += ρ_tmp[ ir, inθ,  iz] * wzr[inz] * wrr[inr] * wθl[inθ]
                            ρ_cell += ρ_tmp[ ir, inθ, inz] * wzl[inz] * wrr[inr] * wθl[inθ]

                            ρ_cell += ρ_tmp[inr,  iθ,  iz] * wzr[inz] * wrl[inr] * wθr[inθ]
                            ρ_cell += ρ_tmp[inr,  iθ, inz] * wzl[inz] * wrl[inr] * wθr[inθ]
                            ρ_cell += ρ_tmp[inr, inθ,  iz] * wzr[inz] * wrl[inr] * wθl[inθ]
                            ρ_cell += ρ_tmp[inr, inθ, inz] * wzl[inz] * wrl[inr] * wθl[inθ]
                        else
                            ρ_cell += ρ_tmp[ ir,  iθ,  iz] * wzr[inz] * 0.5 #wθr[inθ]
                            ρ_cell += ρ_tmp[ ir,  iθ, inz] * wzl[inz] * 0.5 #wθr[inθ]
                            ρ_cell += ρ_tmp[ ir, inθ,  iz] * wzr[inz] * 0.5 #wθl[inθ]
                            ρ_cell += ρ_tmp[ ir, inθ, inz] * wzl[inz] * 0.5 #wθl[inθ]
                        end
                    end
                    if inr > 1
                        wrr_eps::T = ϵ[  ir,  iθ, inz + 1] * wθr[inθ] * wzr[inz]
                        wrr_eps   += ϵ[  ir, inθ, inz + 1] * wθl[inθ] * wzr[inz]
                        wrr_eps   += ϵ[  ir,  iθ, inz ]    * wθr[inθ] * wzl[inz]
                        wrr_eps   += ϵ[  ir, inθ, inz ]    * wθl[inθ] * wzl[inz]
                        # # left weight in r: wrr
                        wrl_eps::T = ϵ[ inr,  iθ, inz + 1] * wθr[inθ] * wzr[inz]
                        wrl_eps   += ϵ[ inr, inθ, inz + 1] * wθl[inθ] * wzr[inz]
                        wrl_eps   += ϵ[ inr,  iθ, inz ]    * wθr[inθ] * wzl[inz]
                        wrl_eps   += ϵ[ inr, inθ, inz ]    * wθl[inθ] * wzl[inz]
                        # right weight in θ: wθr
                        wθr_eps::T = ϵ[ inr,  iθ, inz + 1] * wrl[inr] * wzr[inz]
                        wθr_eps   += ϵ[  ir,  iθ, inz + 1] * wrr[inr] * wzr[inz]
                        wθr_eps   += ϵ[ inr,  iθ, inz ]    * wrl[inr] * wzl[inz]
                        wθr_eps   += ϵ[  ir,  iθ, inz ]    * wrr[inr] * wzl[inz]
                        # left weight in θ: wθl
                        wθl_eps::T = ϵ[ inr, inθ, inz + 1] * wrl[inr] * wzr[inz]
                        wθl_eps   += ϵ[  ir, inθ, inz + 1] * wrr[inr] * wzr[inz]
                        wθl_eps   += ϵ[ inr, inθ, inz ]    * wrl[inr] * wzl[inz]
                        wθl_eps   += ϵ[  ir, inθ, inz ]    * wrr[inr] * wzl[inz]
                        # right weight in z: wzr
                        wzr_eps::T = ϵ[  ir,  iθ, inz + 1] * wrr[inr] * wθr[inθ]
                        wzr_eps   += ϵ[  ir, inθ, inz + 1] * wrr[inr] * wθl[inθ]
                        wzr_eps   += ϵ[ inr,  iθ, inz + 1] * wrl[inr] * wθr[inθ]
                        wzr_eps   += ϵ[ inr, inθ, inz + 1] * wrl[inr] * wθl[inθ]
                        # left weight in z: wzr
                        wzl_eps::T = ϵ[  ir,  iθ,    inz ] * wrr[inr] * wθr[inθ]
                        wzl_eps   += ϵ[  ir, inθ,    inz ] * wrr[inr] * wθl[inθ]
                        wzl_eps   += ϵ[ inr,  iθ,    inz ] * wrl[inr] * wθr[inθ]
                        wzl_eps   += ϵ[ inr, inθ,    inz ] * wrl[inr] * wθl[inθ]

                        volume_weight::T = wrr_eps * Δr_ext_inv[ ir] * mpr[ ir] * Δmpθ[inθ] * Δmpz[inz]
                        volume_weight += wrl_eps * Δr_ext_inv[inr] * mpr[inr] * Δmpθ[inθ] * Δmpz[inz]
                        volume_weight += wθr_eps * r_inv[inr] * Δθ_ext_inv[ iθ] * Δmpr[inr] * Δmpz[inz]
                        volume_weight += wθl_eps * r_inv[inr] * Δθ_ext_inv[inθ] * Δmpr[inr] * Δmpz[inz]
                        volume_weight += wzr_eps * Δz_ext_inv[ iz] * Δmpθ[inθ] * Δmpr_squared[inr]
                        volume_weight += wzl_eps * Δz_ext_inv[inz] * Δmpθ[inθ] * Δmpr_squared[inr]

                        volume_weights[ irbz, iθ, ir, rbi ] = inv(volume_weight)
                        
                        dV::T = Δmpz[inz] * Δmpθ[inθ] * Δmpr_squared[inr]
                        ρ[ irbz, iθ, ir, rbi ] = dV * ρ_cell
                    else
                        wrr_eps = ϵ[  ir,  iθ, inz + 1] * 0.5 * wzr[inz]
                        wrr_eps   += ϵ[  ir, inθ, inz + 1] * 0.5 * wzr[inz]
                        wrr_eps   += ϵ[  ir,  iθ, inz ]    * 0.5 * wzl[inz]
                        wrr_eps   += ϵ[  ir, inθ, inz ]    * 0.5 * wzl[inz]
                        # # left weight in r: wrr
                        wrl_eps = ϵ[ inr,  iθ, inz + 1] * 0.5 * wzr[inz]
                        wrl_eps   += ϵ[ inr, inθ, inz + 1] * 0.5 * wzr[inz]
                        wrl_eps   += ϵ[ inr,  iθ, inz ]    * 0.5 * wzl[inz]
                        wrl_eps   += ϵ[ inr, inθ, inz ]    * 0.5 * wzl[inz]
                        # right weight in θ: wθr
                        wθr_eps = ϵ[ inr,  iθ, inz + 1] * wrl[inr] * wzr[inz]
                        wθr_eps   += ϵ[  ir,  iθ, inz + 1] * wrr[inr] * wzr[inz]
                        wθr_eps   += ϵ[ inr,  iθ, inz ]    * wrl[inr] * wzl[inz]
                        wθr_eps   += ϵ[  ir,  iθ, inz ]    * wrr[inr] * wzl[inz]
                        # left weight in θ: wθl
                        wθl_eps = ϵ[ inr, inθ, inz + 1] * wrl[inr] * wzr[inz]
                        wθl_eps   += ϵ[  ir, inθ, inz + 1] * wrr[inr] * wzr[inz]
                        wθl_eps   += ϵ[ inr, inθ, inz ]    * wrl[inr] * wzl[inz]
                        wθl_eps   += ϵ[  ir, inθ, inz ]    * wrr[inr] * wzl[inz]
                        # right weight in z: wzr
                        wzr_eps = ϵ[  ir,  iθ, inz + 1] * wrr[inr] * 0.5
                        wzr_eps   += ϵ[  ir, inθ, inz + 1] * wrr[inr] * 0.5
                        wzr_eps   += ϵ[ inr,  iθ, inz + 1] * wrl[inr] * 0.5
                        wzr_eps   += ϵ[ inr, inθ, inz + 1] * wrl[inr] * 0.5
                        # left weight in z: wzr
                        wzl_eps = ϵ[  ir,  iθ,    inz ] * wrr[inr] * 0.5
                        wzl_eps   += ϵ[  ir, inθ,    inz ] * wrr[inr] * 0.5
                        wzl_eps   += ϵ[ inr,  iθ,    inz ] * wrl[inr] * 0.5
                        wzl_eps   += ϵ[ inr, inθ,    inz ] * wrl[inr] * 0.5

                        # if inr == 1
                        #     pwwθr = 0.5f0
                        #     pwwθl = 0.5f0
                        #     pwΔmpθ = 2π
                        #     Δθ_ext_inv_r = 0.15915494f0
                        #     Δθ_ext_inv_l = 0.15915494f0
                        # end

                        volume_weight = wrr_eps * Δr_ext_inv[ ir] * mpr[ ir] * 2π * Δmpz[inz]
                        volume_weight += wrl_eps * Δr_ext_inv[inr] * mpr[inr] * Δmpθ[inθ] * Δmpz[inz]
                        volume_weight += wθr_eps * r_inv[inr] * 0.15915494f0 * Δmpr[inr] * Δmpz[inz]
                        volume_weight += wθl_eps * r_inv[inr] * 0.15915494f0 * Δmpr[inr] * Δmpz[inz]
                        volume_weight += wzr_eps * Δz_ext_inv[ iz] * 2π * Δmpr_squared[inr]
                        volume_weight += wzl_eps * Δz_ext_inv[inz] * 2π * Δmpr_squared[inr]

                        volume_weights[ irbz, iθ, ir, rbi ] = inv(volume_weight)
                        
                        dV = Δmpz[inz] * 2π * Δmpr_squared[inr]
                        ρ[ irbz, iθ, ir, rbi ] = dV * ρ_cell
                    end
                end
            end
        end
        
        potential::Array{T, 3} = ismissing(potential_array) ? zeros(T, size(grid)...) : potential_array
        pointtypes::Array{PointType, 3} = ones(PointType, size(grid)...)
        set_pointtypes_and_fixed_potentials!( pointtypes, potential, grid, ssd, weighting_potential_channel_idx = weighting_potential_channel_idx  )
        rbpotential::Array{T, 4}  = RBExtBy2Array( potential, grid )
        rbpointtypes::Array{T, 4} = RBExtBy2Array( pointtypes, grid )
        potential = clear(potential)
        pointtypes = clear(pointtypes)
    end # @inbounds
        
    fssrb::PotentialSimulationSetupRB{T, 3, 4, :Cylindrical} = PotentialSimulationSetupRB{T, 3, 4, :Cylindrical}(  
        grid, 
        rbpotential,
        rbpointtypes,
        volume_weights,
        ρ,
        ϵ,
        geom_weights,
        sor_const,
        bias_voltage,
        maximum_applied_potential,
        minimum_applied_potential,
        depletion_handling_potential_limit,
        bulk_is_ptype,
        grid_boundary_factors
     )
    return fssrb
end

function Grid(fssrb::PotentialSimulationSetupRB{T, N1, N2, S})::Grid{T, N1, S} where {T, N1, N2, S}
    return fssrb.grid
end

function ElectricPotentialArray(fssrb::PotentialSimulationSetupRB{T, 3, 4, :Cylindrical})::Array{T, 3} where {T}
    pot::Array{T, 3} = Array{T, 3}(undef, size(fssrb.grid))
    for iz in axes(pot, 3)
        irbz::Int = rbidx(iz)
        for iθ in axes(pot, 2)
            irbθ::Int = iθ + 1
            idxsum::Int = iz + iθ
            for ir in axes(pot, 1)
                irbr::Int = ir + 1
                rbi::Int = iseven(idxsum + ir) ? rb_even::Int : rb_odd::Int
                pot[ir, iθ, iz] = fssrb.potential[ irbz, irbθ, irbr, rbi ]
            end
        end
    end
    return pot
end


function PointTypeArray(fssrb::PotentialSimulationSetupRB{T, 3, 4, :Cylindrical})::Array{PointType, 3} where {T}
    pointtypes::Array{PointType, 3} = zeros(PointType, size(fssrb.grid))
    for iz in axes(pointtypes, 3)
        irbz::Int = rbidx(iz)
        for iθ in axes(pointtypes, 2)
            irbθ::Int = iθ + 1
            idxsum::Int = iz + iθ
            for ir in axes(pointtypes, 1)
                irbr::Int = ir + 1
                rbi::Int = iseven(idxsum + ir) ? rb_even::Int : rb_odd::Int
                pointtypes[ir, iθ, iz] = fssrb.pointtypes[irbz, irbθ, irbr, rbi ]
            end
        end
    end
    return pointtypes
end


function PointTypes(fss::PotentialSimulationSetup{T, N, S})::PointTypes{T, N, S} where {T, N, S}
    return PointTypes{T, N, S}( fss.pointtypes, fss.grid )
end


function ChargeDensityArray(fssrb::PotentialSimulationSetupRB{T, 3, 4, :Cylindrical})::Array{T} where {T}
    ρ::Array{T, 3} = zeros(T, size(fssrb.grid))
    for iz in axes(ρ, 3)
        irbz::Int = rbidx(iz)
        Δmpz::T = fssrb.geom_weights[3].weights[3, iz] 
        for iθ in axes(ρ, 2)
            irbθ::Int = iθ + 1
            idxsum::Int = iz + iθ
            Δmpθ::T = fssrb.geom_weights[2].weights[3, iθ]
            Δmpzθ::T = Δmpz * Δmpθ
            for ir in axes(ρ, 1)
                irbr::Int = ir + 1
                rbi::Int = iseven(idxsum + ir) ? rb_even::Int : rb_odd::Int
                dV::T = fssrb.geom_weights[1].weights[6, ir] * Δmpzθ  #Δmpz[inz] * Δmpθ[inθ] * Δmpr_squared[inr]
                if ir == 1
                    dV = dV * 2π / Δmpθ
                end
                ρ[ir, iθ, iz] = fssrb.ρ[irbz, irbθ, irbr, rbi ] / dV
            end
        end
    end
    return ρ
end


function DielektrikumDistributionArray(fssrb::PotentialSimulationSetupRB{T, 3, 4, :Cylindrical})::Array{T, 3} where {T}
    return fssrb.ϵ 
end

