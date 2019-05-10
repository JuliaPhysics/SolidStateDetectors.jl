function PotentialSimulationSetupRB(ssd::SolidStateDetector{T, :cartesian}, grid::Grid{T, 3, :cartesian} = Grid(ssd),
                potential_array::Union{Missing, Array{T, 3}} = missing; weighting_potential_contact_id::Union{Missing, Int} = missing) where {T}#::PotentialSimulationSetupRB{T} where {T}
    is_weighting_potential::Bool = !ismissing(weighting_potential_contact_id)

    @inbounds begin
        begin # Geometrical weights of the Axes
            nx, ny, nz = size(grid)

            # X-axis
            ax_x::Vector{T} = collect(grid.axes[1]) # real grid points/ticks -> the potential at these ticks are going to be calculated
            Δx::Vector{T} = diff(ax_x) # difference between real ticks
            x_ext::Vector{T} = get_extended_ticks(grid.axes[1])
            Δx_ext::Vector{T} = diff(x_ext)
            Δx_ext_inv::Vector{T} = inv.(Δx_ext)
            mpx::Vector{T} = midpoints(x_ext)
            mpx_inv::Vector{T} = inv.(mpx)
            Δmpx::Vector{T} = diff(mpx)
            Δmpx_inv::Vector{T} = inv.(Δmpx)
            Δhmpxr::Vector{T} = mpx[2:end] - ax_x # distances between midpoints and real grid points (half distances -> h), needed for weights wrr, wrl, ... & dV (volume element)
            Δhmpxl::Vector{T} = ax_x - mpx[1:end - 1]
            wxr::Vector{T} = Δmpx_inv .* Δhmpxr # weights for epislon_r adding
            wxl::Vector{T} = Δmpx_inv .* Δhmpxl
            wx::Array{T, 2} = zeros(T, 4, length(wxr) + 1)
            wx[1, 1:length(wxr)] = wxr
            wx[2, 1:length(wxr)] = wxl
            wx[3, 1:length(wxr)] = Δmpx
            wx[4, :] = Δx_ext_inv
            gw_x::GeometricalCartesianAxisWeights{T} = GeometricalCartesianAxisWeights{T}( wx ) # Weights needed for Field Simulation loop
            x_ext_mid::T = (x_ext[end - 1] - x_ext[2]) / 2
            grid_boundary_factor_x_right::T = abs((x_ext[end - 1] - x_ext_mid) / (x_ext[end] - x_ext_mid))
            grid_boundary_factor_x_left::T = abs((x_ext[2] - x_ext_mid) / (x_ext[1] - x_ext_mid))

            # Y-axis
            ax_y::Vector{T} = collect(grid.axes[2]) # real grid points/ticks -> the potential at these ticks are going to be calculated
            Δy::Vector{T} = diff(ax_y) # difference between real ticks
            y_ext::Vector{T} = get_extended_ticks(grid.axes[2])
            Δy_ext::Vector{T} = diff(y_ext)
            Δy_ext_inv::Vector{T} = inv.(Δy_ext)
            mpy::Vector{T} = midpoints(y_ext)
            mpy_inv::Vector{T} = inv.(mpy)
            Δmpy::Vector{T} = diff(mpy)
            Δmpy_inv::Vector{T} = inv.(Δmpy)
            Δhmpyr::Vector{T} = mpy[2:end] - ax_y # distances between midpoints and real grid points (half distances -> h), needed for weights wrr, wrl, ... & dV (volume element)
            Δhmpyl::Vector{T} = ax_y - mpy[1:end - 1]
            wyr::Vector{T} = Δmpy_inv .* Δhmpyr # weights for epislon_r adding
            wyl::Vector{T} = Δmpy_inv .* Δhmpyl
            wy::Array{T, 2} = zeros(T, 4, length(wyr) + 1)
            wy[1, 1:length(wyr)] = wyr
            wy[2, 1:length(wyr)] = wyl
            wy[3, 1:length(wyr)] = Δmpy
            wy[4, :] = Δy_ext_inv
            gw_y::GeometricalCartesianAxisWeights{T} = GeometricalCartesianAxisWeights{T}( wy ) # Weights needed for Field Simulation loop
            y_ext_mid::T = (y_ext[end - 1] - y_ext[2]) / 2
            grid_boundary_factor_y_right::T = abs((y_ext[end - 1] - y_ext_mid) / (y_ext[end] - y_ext_mid))
            grid_boundary_factor_y_left::T = abs((y_ext[2] - y_ext_mid) / (y_ext[1] - y_ext_mid))

            # Z-axis
            ax_z::Vector{T} = collect(grid.axes[3]) # real grid points/ticks -> the potential at these ticks are going to be calculated
            Δz::Vector{T} = diff(ax_z) # difference between real ticks
            z_ext::Vector{T} = get_extended_ticks(grid.axes[3])
            Δz_ext::Vector{T} = diff(z_ext)
            Δz_ext_inv::Vector{T} = inv.(Δz_ext)
            mpz::Vector{T} = midpoints(z_ext)
            mpz_inv::Vector{T} = inv.(mpz)
            Δmpz::Vector{T} = diff(mpz)
            Δmpz_inv::Vector{T} = inv.(Δmpz)
            Δhmpzr::Vector{T} = mpz[2:end] - ax_z # distances between midpoints and real grid points (half distances -> h), needed for weights wrr, wrl, ... & dV (volume element)
            Δhmpzl::Vector{T} = ax_z - mpz[1:end - 1]
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

            geom_weights::NTuple{3, AbstractGeometricalAxisWeights{T}} = (gw_x, gw_y, gw_z) # Weights needed for Field Simulation loop
            grid_boundary_factors::NTuple{3, NTuple{2, T}} = (  (grid_boundary_factor_x_left, grid_boundary_factor_x_right),
                                                                (grid_boundary_factor_y_left, grid_boundary_factor_y_right),
                                                                (grid_boundary_factor_z_left, grid_boundary_factor_z_right))
        end

        # detector_material_ϵ_r::T = ssd.material_detector.ϵ_r
        # environment_material_ϵ_r::T = ssd.material_environment.ϵ_r

        bulk_is_ptype::Bool = if length(ssd.semiconductors) > 0
            ssd.semiconductors[1].bulk_type == :ptype ? true : false
        else
            true
        end
        contact_bias_voltages::Vector{T} = if length(ssd.contacts) > 0 
            T[contact.potential for contact in ssd.contacts]
        else
            T[0]
        end
        minimum_applied_potential::T = minimum(contact_bias_voltages)
        maximum_applied_potential::T = maximum(contact_bias_voltages)
        bias_voltage::T = maximum_applied_potential - minimum_applied_potential
        depletion_handling_potential_limit::T = -bias_voltage
        sor_consts::Vector{T} = [1]

        ϵ = Array{T, 3}(undef, length(mpx), length(mpy), length(mpz))
        ρ_tmp = Array{T, 3}(undef, length(mpx), length(mpy), length(mpz))
        for iz in 1:size(ϵ, 3)
            pos_z::T = mpz[iz]
            for iy in 1:size(ϵ, 2)
                pos_y::T = mpy[iy]
                for ix in 1:size(ϵ, 1)
                    pos_x = mpx[ix]
                    pt::CartesianPoint{T} = CartesianPoint{T}(pos_x, pos_y, pos_z)
                    ρ_tmp[ix, iy, iz]::T, ϵ[ix, iy, iz]::T = get_ρ_and_ϵ(pt, ssd)
                end
            end
        end
        ϵ0_inv::T = inv(ϵ0)
        ρ_tmp *= ϵ0_inv

        volume_weights::Array{T, 4} = RBExtBy2Array(T, grid)
        ρ::Array{T, 4} = RBExtBy2Array(T, grid)
        for iz in range(2, stop = length(z_ext) - 1)
            inz::Int = iz - 1
            pos_z::T = z_ext[iz]
            for iy in range(2, stop = length(y_ext) - 1)
                iny::Int = iy - 1
                pos_y::T = y_ext[iy]
                for ix in range(2, stop = length(x_ext) - 1)
                    inx::Int = ix - 1;
                    irbx::Int = rbidx(inx)
                    pos_x::T = x_ext[ix]

                    rbi::Int = iseven(inx + iny + inz) ? rb_even::Int : rb_odd::Int

                    ρ_cell::T = 0
                    if !is_weighting_potential
                        ρ_cell += ρ_tmp[ ix,  iy,  iz] * wxr[inx] * wyr[iny] * wzr[inz]
                        ρ_cell += ρ_tmp[ ix,  iy, inz] * wxr[inx] * wyr[iny] * wzl[inz]
                        ρ_cell += ρ_tmp[ ix, iny,  iz] * wxr[inx] * wyl[iny] * wzr[inz]
                        ρ_cell += ρ_tmp[ ix, iny, inz] * wxr[inx] * wyl[iny] * wzl[inz]

                        ρ_cell += ρ_tmp[inx,  iy,  iz] * wxl[inx] * wyr[iny] * wzr[inz]
                        ρ_cell += ρ_tmp[inx,  iy, inz] * wxl[inx] * wyr[iny] * wzl[inz]
                        ρ_cell += ρ_tmp[inx, iny,  iz] * wxl[inx] * wyl[iny] * wzr[inz]
                        ρ_cell += ρ_tmp[inx, iny, inz] * wxl[inx] * wyl[iny] * wzl[inz]
                    end

                    # right weight in x: wxr
                    wxr_eps::T = ϵ[  ix,  iy, inz + 1] * wyr[iny] * wzr[inz]
                    wxr_eps   += ϵ[  ix, iny, inz + 1] * wyl[iny] * wzr[inz]
                    wxr_eps   += ϵ[  ix,  iy, inz ]    * wyr[iny] * wzl[inz]
                    wxr_eps   += ϵ[  ix, iny, inz ]    * wyl[iny] * wzl[inz]
                    # left weight in x: wxr
                    wxl_eps::T = ϵ[ inx,  iy, inz + 1] * wyr[iny] * wzr[inz]
                    wxl_eps   += ϵ[ inx, iny, inz + 1] * wyl[iny] * wzr[inz]
                    wxl_eps   += ϵ[ inx,  iy, inz ]    * wyr[iny] * wzl[inz]
                    wxl_eps   += ϵ[ inx, iny, inz ]    * wyl[iny] * wzl[inz]
                    # right weight in y: wyr
                    wyr_eps::T = ϵ[ inx,  iy, inz + 1] * wxl[inx] * wzr[inz]
                    wyr_eps   += ϵ[  ix,  iy, inz + 1] * wxr[inx] * wzr[inz]
                    wyr_eps   += ϵ[ inx,  iy, inz ]    * wxl[inx] * wzl[inz]
                    wyr_eps   += ϵ[  ix,  iy, inz ]    * wxr[inx] * wzl[inz]
                    # left weight in y: wyl
                    wyl_eps::T = ϵ[ inx, iny, inz + 1] * wxl[inx] * wzr[inz]
                    wyl_eps   += ϵ[  ix, iny, inz + 1] * wxr[inx] * wzr[inz]
                    wyl_eps   += ϵ[ inx, iny, inz ]    * wxl[inx] * wzl[inz]
                    wyl_eps   += ϵ[  ix, iny, inz ]    * wxr[inx] * wzl[inz]
                    # right weight in z: wzr
                    wzr_eps::T = ϵ[  ix,  iy, inz + 1] * wxr[inx] * wyr[iny]
                    wzr_eps   += ϵ[  ix, iny, inz + 1] * wxr[inx] * wyl[iny]
                    wzr_eps   += ϵ[ inx,  iy, inz + 1] * wxl[inx] * wyr[iny]
                    wzr_eps   += ϵ[ inx, iny, inz + 1] * wxl[inx] * wyl[iny]
                    # left weight in z: wzr
                    wzl_eps::T = ϵ[  ix,  iy,    inz ] * wxr[inx] * wyr[iny]
                    wzl_eps   += ϵ[  ix, iny,    inz ] * wxr[inx] * wyl[iny]
                    wzl_eps   += ϵ[ inx,  iy,    inz ] * wxl[inx] * wyr[iny]
                    wzl_eps   += ϵ[ inx, iny,    inz ] * wxl[inx] * wyl[iny]

                    volume_weight::T = wxr_eps * Δx_ext_inv[ ix] * Δmpy[iny] * Δmpz[inz]
                    volume_weight   += wxl_eps * Δx_ext_inv[inx] * Δmpy[iny] * Δmpz[inz]
                    volume_weight   += wyr_eps * Δy_ext_inv[ iy] * Δmpx[inx] * Δmpz[inz]
                    volume_weight   += wyl_eps * Δy_ext_inv[iny] * Δmpx[inx] * Δmpz[inz]
                    volume_weight   += wzr_eps * Δz_ext_inv[ iz] * Δmpx[inx] * Δmpy[iny]
                    volume_weight   += wzl_eps * Δz_ext_inv[inz] * Δmpx[inx] * Δmpy[iny]

                    volume_weights[ irbx, iy, iz, rbi ] = inv(volume_weight)

                    dV::T = Δmpx[inx] * Δmpy[iny] * Δmpz[inz]
                    ρ[ irbx, iy, iz, rbi ] = dV * ρ_cell
                end
            end
        end

        potential::Array{T, 3} = ismissing(potential_array) ? zeros(T, size(grid)...) : potential_array
        pointtypes::Array{PointType, 3} = ones(PointType, size(grid)...)
        set_pointtypes_and_fixed_potentials!( pointtypes, potential, grid, ssd, weighting_potential_contact_id = weighting_potential_contact_id  )
        rbpotential::Array{T, 4}  = RBExtBy2Array( potential, grid )
        rbpointtypes::Array{T, 4} = RBExtBy2Array( pointtypes, grid )
        potential = clear(potential)
        pointtypes = clear(pointtypes)
    end # @inbounds

    fssrb::PotentialSimulationSetupRB{T, 3, 4, :cartesian} = PotentialSimulationSetupRB{T, 3, 4, :cartesian}(
        grid,
        rbpotential,
        rbpointtypes,
        volume_weights,
        ρ,
        ϵ,
        geom_weights,
        sor_consts,
        bias_voltage,
        maximum_applied_potential,
        minimum_applied_potential,
        depletion_handling_potential_limit,
        bulk_is_ptype,
        grid_boundary_factors
     )
    return fssrb
end


function ElectricPotentialArray(fssrb::PotentialSimulationSetupRB{T, 3, 4, :cartesian})::Array{T, 3} where {T}
    pot::Array{T, 3} = Array{T, 3}(undef, size(fssrb.grid))
    for iz in axes(pot, 3)
        irbz::Int = iz + 1
        for iy in axes(pot, 2)
            irby::Int = iy + 1
            idxsum::Int = iz + iy
            for ix in axes(pot, 1)
                irbx::Int = rbidx(ix)
                rbi::Int = iseven(idxsum + ix) ? rb_even::Int : rb_odd::Int
                pot[ix, iy, iz] = fssrb.potential[ irbx, irby, irbz, rbi ]
            end
        end
    end
    return pot
end


function PointTypeArray(fssrb::PotentialSimulationSetupRB{T, 3, 4, :cartesian})::Array{PointType, 3} where {T}
    pointtypes::Array{PointType, 3} = zeros(PointType, size(fssrb.grid))
    for iz in axes(pointtypes, 3)
        irbz::Int = iz + 1
        for iy in axes(pointtypes, 2)
            irby::Int = iy + 1
            idxsum::Int = iz + iy
            for ix in axes(pointtypes, 1)
                irbx::Int = rbidx(ix)
                rbi::Int = iseven(idxsum + ix) ? rb_even::Int : rb_odd::Int
                pointtypes[ix, iy, iz] = fssrb.pointtypes[irbx, irby, irbz, rbi ]
            end
        end
    end
    return pointtypes
end



function ChargeDensityArray(fssrb::PotentialSimulationSetupRB{T, 3, 4, :cartesian})::Array{T} where {T}
    ρ::Array{T, 3} = zeros(T, size(fssrb.grid))
    for iz in axes(ρ, 3)
        irbz::Int = iz + 1
        Δmpz::T = fssrb.geom_weights[3].weights[3, iz]
        for iy in axes(ρ, 2)
            irby::Int = iy + 1
            idxsum::Int = iz + iy
            Δmpy::T = fssrb.geom_weights[2].weights[3, iy]
            Δmpzy::T = Δmpz * Δmpy
            for ix in axes(ρ, 1)
                irbx::Int = rbidx(ix)
                rbi::Int = iseven(idxsum + ix) ? rb_even::Int : rb_odd::Int
                dV::T =  Δmpzy * fssrb.geom_weights[1].weights[3, ix]

                ρ[ix, iy, iz] = fssrb.ρ[irbx, irby, irbz, rbi ] / dV
            end
        end
    end
    return ρ
end
