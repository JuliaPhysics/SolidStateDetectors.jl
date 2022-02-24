function set_passive_or_contact_points(point_types::Array{PointType, 3}, potential::Array{T, 3},
                        grid::CartesianGrid3D{T}, obj, pot::T, use_nthreads::Int = 1) where {T}
    if !isnan(pot)
        @onthreads 1:use_nthreads for iz in workpart(axes(potential, 3), 1:use_nthreads, Base.Threads.threadid())
            @inbounds for iy in axes(potential, 2)
                for ix in axes(potential, 1)
                    pt::CartesianPoint{T} = CartesianPoint{T}( grid.axes[1].ticks[ix], grid.axes[2].ticks[iy], grid.axes[3].ticks[iz] )
                    if pt in obj
                        potential[ ix, iy, iz ] = pot
                        point_types[ ix, iy, iz ] = zero(PointType)
                    end
                end
            end
        end
    end
    nothing
end

function set_point_types_and_fixed_potentials!(point_types::Array{PointType, 3}, potential::Array{T, 3},
        grid::CartesianGrid3D{T}, det::SolidStateDetector{T}; 
        weighting_potential_contact_id::Union{Missing, Int} = missing,
        use_nthreads::Int = Base.Threads.nthreads(),
        not_only_paint_contacts::Val{NotOnlyPaintContacts} = Val{true}(),
        paint_contacts::Val{PaintContacts} = Val{true}())::Nothing where {T <: SSDFloat, NotOnlyPaintContacts, PaintContacts}

    @onthreads 1:use_nthreads for iz in workpart(axes(potential, 3), 1:use_nthreads, Base.Threads.threadid())
        @inbounds for iy in axes(potential, 2)
            for ix in axes(potential, 1)
                pt::CartesianPoint{T} = CartesianPoint{T}( grid.axes[1].ticks[ix], grid.axes[2].ticks[iy], grid.axes[3].ticks[iz] )
                if in(pt, det.semiconductor)
                    point_types[ ix, iy, iz ] += pn_junction_bit
                end
            end
        end
    end   
    if !ismissing(det.passives)
        for passive in det.passives
            pot::T = ismissing(weighting_potential_contact_id) ? passive.potential : zero(T)
            set_passive_or_contact_points(point_types, potential, grid, passive.geometry, pot, use_nthreads)                              
        end
    end
    if NotOnlyPaintContacts
        for contact in det.contacts
            pot::T = ismissing(weighting_potential_contact_id) ? contact.potential : contact.id == weighting_potential_contact_id
            set_passive_or_contact_points(point_types, potential, grid, contact.geometry, pot, use_nthreads)
        end
    end
    if PaintContacts
        for contact in det.contacts
            pot::T = ismissing(weighting_potential_contact_id) ? contact.potential : contact.id == weighting_potential_contact_id
            fs = ConstructiveSolidGeometry.surfaces(contact.geometry)
            for face in fs
                paint!(point_types, potential, face, contact.geometry, pot, grid)
            end
        end
    end
    nothing
end

function fill_ρ_and_ϵ!(ϵ::Array{T}, ρ_tmp::Array{T}, q_eff_fix_tmp::Array{T}, 
    ::Type{Cartesian}, mpz::Vector{T}, mpy::Vector{T}, mpx::Vector{T}, use_nthreads::Int, obj) where {T}
    @inbounds begin
        @onthreads 1:use_nthreads for iz in workpart(axes(ϵ, 3), 1:use_nthreads, Base.Threads.threadid())
            for iy in axes(ϵ, 2)
                for ix in axes(ϵ, 1)
                    pt::CartesianPoint{T} = CartesianPoint{T}(mpx[ix], mpy[iy], mpz[iz])
                    if pt in obj
                        ρ_tmp[ix, iy, iz]::T, ϵ[ix, iy, iz]::T, q_eff_fix_tmp[ix, iy, iz]::T = get_ρ_and_ϵ(pt, obj)
                    end
                end
            end
        end
    end
    nothing
end

function PotentialCalculationSetup(det::SolidStateDetector{T}, grid::CartesianGrid3D{T}, 
                medium::NamedTuple = material_properties[materials["vacuum"]],
                potential_array::Union{Missing, Array{T, 3}} = missing,
                imp_scale::Union{Missing, ElectricPotential{T, 3}} = missing;
                weighting_potential_contact_id::Union{Missing, Int} = missing,
                point_types = missing,
                use_nthreads::Int = Base.Threads.nthreads(),
                sor_consts::T = T(1), not_only_paint_contacts::Bool = true, paint_contacts::Bool = true ) where {T}
    
    is_weighting_potential::Bool = !ismissing(weighting_potential_contact_id)
    depletion_handling::Bool = is_weighting_potential && !ismissing(point_types)
    
    @inbounds begin
        begin # Geometrical weights of the Axes
            # X-axis
            ax_x::Vector{T} = collect(grid.axes[1]) # real grid points/ticks -> the potential at these ticks are going to be calculated
            x_ext::Vector{T} = get_extended_ticks(grid.axes[1])
            Δx_ext::Vector{T} = diff(x_ext)
            Δx_ext_inv::Vector{T} = inv.(Δx_ext)
            mpx::Vector{T} = midpoints(x_ext)
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
            y_ext::Vector{T} = get_extended_ticks(grid.axes[2])
            Δy_ext::Vector{T} = diff(y_ext)
            Δy_ext_inv::Vector{T} = inv.(Δy_ext)
            mpy::Vector{T} = midpoints(y_ext)
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
            z_ext::Vector{T} = get_extended_ticks(grid.axes[3])
            Δz_ext::Vector{T} = diff(z_ext)
            Δz_ext_inv::Vector{T} = inv.(Δz_ext)
            mpz::Vector{T} = midpoints(z_ext)
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

            geom_weights::NTuple{3, AbstractGeometricalAxisWeights{T}} = (gw_z, gw_y, gw_x) # Weights needed for Field Simulation loop
            grid_boundary_factors::NTuple{3, NTuple{2, T}} = (  (grid_boundary_factor_x_left, grid_boundary_factor_x_right),
                                                                (grid_boundary_factor_y_left, grid_boundary_factor_y_right),
                                                                (grid_boundary_factor_z_left, grid_boundary_factor_z_right))
        end

        contact_bias_voltages::Vector{T} = if length(det.contacts) > 0 
            T[contact.potential for contact in det.contacts]
        else
            T[0]
        end
        minimum_applied_potential::T = minimum(contact_bias_voltages)
        maximum_applied_potential::T = maximum(contact_bias_voltages)
        bias_voltage::T = maximum_applied_potential - minimum_applied_potential
        sor_consts = [sor_consts]

        medium_ϵ_r::T = medium.ϵ_r
        ϵ = fill(medium_ϵ_r, length(mpx), length(mpy), length(mpz))
        ρ_tmp = zeros(T, length(mpx), length(mpy), length(mpz))
        q_eff_fix_tmp = zeros(T, length(mpx), length(mpy), length(mpz))
        fill_ρ_and_ϵ!(ϵ, ρ_tmp, q_eff_fix_tmp, Cartesian, mpz, mpy, mpx, use_nthreads, det.semiconductor)
        if !ismissing(det.passives)
            for passive in det.passives
                fill_ρ_and_ϵ!(ϵ, ρ_tmp, q_eff_fix_tmp, Cartesian, mpz, mpy, mpx, use_nthreads, passive)
            end
        end
        if depletion_handling
            for iz in axes(ϵ, 3)
                for iy in axes(ϵ, 2)
                    for ix in axes(ϵ, 1)
                        pt::CartesianPoint{T} = CartesianPoint{T}(mpx[ix], mpy[iy], mpz[iz])
                        ig = find_closest_gridpoint(pt, point_types.grid)
                        if is_undepleted_point_type(point_types.data[ig...]) 
                            ϵ[ix, iy, iz] *= scaling_factor_for_permittivity_in_undepleted_region(det.semiconductor) * (1 - imp_scale.data[ig...])
                        end
                    end
                end
            end
        end  

        ϵ0_inv::T = inv(ϵ0)
        ρ_tmp *= ϵ0_inv
        q_eff_fix_tmp *= ϵ0_inv

        volume_weights::Array{T, 4} = RBExtBy2Array(T, grid)
        q_eff_imp::Array{T, 4} = RBExtBy2Array(T, grid)
        q_eff_fix::Array{T, 4} = RBExtBy2Array(T, grid)
        for iz in range(2, stop = length(z_ext) - 1)
            inz::Int = iz - 1
            for iy in range(2, stop = length(y_ext) - 1)
                iny::Int = iy - 1
                for ix in range(2, stop = length(x_ext) - 1)
                    inx::Int = ix - 1;
                    irbx::Int = rbidx(inx)

                    rbi::Int = iseven(inx + iny + inz) ? rb_even::Int : rb_odd::Int

                    ρ_cell::T = 0
                    q_eff_fix_cell::T = 0
                    if !is_weighting_potential
                        ρ_cell += ρ_tmp[ ix,  iy,  iz] * wxr[inx] * wyr[iny] * wzr[inz]
                        ρ_cell += ρ_tmp[ ix,  iy, inz] * wxr[inx] * wyr[iny] * wzl[inz]
                        ρ_cell += ρ_tmp[ ix, iny,  iz] * wxr[inx] * wyl[iny] * wzr[inz]
                        ρ_cell += ρ_tmp[ ix, iny, inz] * wxr[inx] * wyl[iny] * wzl[inz]

                        ρ_cell += ρ_tmp[inx,  iy,  iz] * wxl[inx] * wyr[iny] * wzr[inz]
                        ρ_cell += ρ_tmp[inx,  iy, inz] * wxl[inx] * wyr[iny] * wzl[inz]
                        ρ_cell += ρ_tmp[inx, iny,  iz] * wxl[inx] * wyl[iny] * wzr[inz]
                        ρ_cell += ρ_tmp[inx, iny, inz] * wxl[inx] * wyl[iny] * wzl[inz]

                        q_eff_fix_cell += q_eff_fix_tmp[ ix,  iy,  iz] * wxr[inx] * wyr[iny] * wzr[inz]
                        q_eff_fix_cell += q_eff_fix_tmp[ ix,  iy, inz] * wxr[inx] * wyr[iny] * wzl[inz]
                        q_eff_fix_cell += q_eff_fix_tmp[ ix, iny,  iz] * wxr[inx] * wyl[iny] * wzr[inz]
                        q_eff_fix_cell += q_eff_fix_tmp[ ix, iny, inz] * wxr[inx] * wyl[iny] * wzl[inz]

                        q_eff_fix_cell += q_eff_fix_tmp[inx,  iy,  iz] * wxl[inx] * wyr[iny] * wzr[inz]
                        q_eff_fix_cell += q_eff_fix_tmp[inx,  iy, inz] * wxl[inx] * wyr[iny] * wzl[inz]
                        q_eff_fix_cell += q_eff_fix_tmp[inx, iny,  iz] * wxl[inx] * wyl[iny] * wzr[inz]
                        q_eff_fix_cell += q_eff_fix_tmp[inx, iny, inz] * wxl[inx] * wyl[iny] * wzl[inz]
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
                    q_eff_imp[ irbx, iy, iz, rbi ] = dV * ρ_cell
                    q_eff_fix[ irbx, iy, iz, rbi ] = dV * q_eff_fix_cell
                end
            end
        end

        potential = ismissing(potential_array) ? zeros(T, size(grid)...) : potential_array
        point_types = ones(PointType, size(grid)...)
        set_point_types_and_fixed_potentials!( point_types, potential, grid, det, 
                weighting_potential_contact_id = weighting_potential_contact_id,
                use_nthreads = use_nthreads,
                not_only_paint_contacts = Val(not_only_paint_contacts), 
                paint_contacts = Val(paint_contacts)  )
        rbpotential = RBExtBy2Array( potential, grid )
        rbpoint_types = RBExtBy2Array( point_types, grid )
        
        imp_scale = ismissing(imp_scale) || is_weighting_potential ? zeros(T, size(q_eff_imp)) : RBExtBy2Array(imp_scale.data, grid)
        for i in eachindex(imp_scale)
            imp_scale[i] = is_pn_junction_point_type(rbpoint_types[i])
        end
    end # @inbounds

    pcs = PotentialCalculationSetup(
        grid,
        rbpotential,
        rbpoint_types,
        volume_weights,
        q_eff_imp,
        imp_scale,
        q_eff_fix,
        ϵ,
        broadcast(gw -> gw.weights, geom_weights),
        sor_consts,
        bias_voltage,
        maximum_applied_potential,
        minimum_applied_potential,
        grid_boundary_factors
    )
    return pcs
end


function ElectricPotentialArray(pcs::PotentialCalculationSetup{T, Cartesian,  3, Array{T, 3}})::Array{T, 3} where {T}
    pot::Array{T, 3} = Array{T, 3}(undef, size(pcs.grid))
    for iz in axes(pot, 3)
        irbz::Int = iz + 1
        for iy in axes(pot, 2)
            irby::Int = iy + 1
            idxsum::Int = iz + iy
            for ix in axes(pot, 1)
                irbx::Int = rbidx(ix)
                rbi::Int = iseven(idxsum + ix) ? rb_even::Int : rb_odd::Int
                pot[ix, iy, iz] = pcs.potential[ irbx, irby, irbz, rbi ]
            end
        end
    end
    return pot
end

function ImpurityScale(pcs::PotentialCalculationSetup{T, Cartesian, 3, Array{T, 3}})::Array{T, 3} where {T}
    s::Array{T, 3} = Array{T, 3}(undef, size(pcs.grid))
    for iz in axes(s, 3)
        irbz::Int = iz + 1
        for iy in axes(s, 2)
            irby::Int = iy + 1
            idxsum::Int = iz + iy
            for ix in axes(s, 1)
                irbx::Int = rbidx(ix)
                rbi::Int = iseven(idxsum + ix) ? rb_even::Int : rb_odd::Int
                s[ix, iy, iz] = pcs.imp_scale[ irbx, irby, irbz, rbi ]
            end
        end
    end
    return s
end

function PointTypeArray(pcs::PotentialCalculationSetup{T, Cartesian,  3, Array{T, 3}})::Array{PointType, 3} where {T}
    point_types::Array{PointType, 3} = zeros(PointType, size(pcs.grid))
    for iz in axes(point_types, 3)
        irbz::Int = iz + 1
        for iy in axes(point_types, 2)
            irby::Int = iy + 1
            idxsum::Int = iz + iy
            for ix in axes(point_types, 1)
                irbx::Int = rbidx(ix)
                rbi::Int = iseven(idxsum + ix) ? rb_even::Int : rb_odd::Int
                point_types[ix, iy, iz] = pcs.point_types[irbx, irby, irbz, rbi ]
            end
        end
    end
    return point_types
end



function EffectiveChargeDensityArray(pcs::PotentialCalculationSetup{T, Cartesian,  3, Array{T, 3}})::Array{T} where {T}
    ρ::Array{T, 3} = zeros(T, size(pcs.grid))
    for iz in axes(ρ, 3)
        irbz::Int = iz + 1
        Δmpz::T = pcs.geom_weights[1][3, iz]
        for iy in axes(ρ, 2)
            irby::Int = iy + 1
            idxsum::Int = iz + iy
            Δmpy::T = pcs.geom_weights[2][3, iy]
            Δmpzy::T = Δmpz * Δmpy
            for ix in axes(ρ, 1)
                irbx::Int = rbidx(ix)
                rbi::Int = iseven(idxsum + ix) ? rb_even::Int : rb_odd::Int
                dV::T =  Δmpzy * pcs.geom_weights[3][3, ix]

                ρ[ix, iy, iz] = pcs.q_eff_imp[irbx, irby, irbz, rbi ] / dV
            end
        end
    end
    return ρ
end
function FixedEffectiveChargeDensityArray(pcs::PotentialCalculationSetup{T, Cartesian,  3, Array{T, 3}})::Array{T} where {T}
    ρ::Array{T, 3} = zeros(T, size(pcs.grid))
    for iz in axes(ρ, 3)
        irbz::Int = iz + 1
        Δmpz::T = pcs.geom_weights[1][3, iz]
        for iy in axes(ρ, 2)
            irby::Int = iy + 1
            idxsum::Int = iz + iy
            Δmpy::T = pcs.geom_weights[2][3, iy]
            Δmpzy::T = Δmpz * Δmpy
            for ix in axes(ρ, 1)
                irbx::Int = rbidx(ix)
                rbi::Int = iseven(idxsum + ix) ? rb_even::Int : rb_odd::Int
                dV::T =  Δmpzy * pcs.geom_weights[3][3, ix]

                ρ[ix, iy, iz] = pcs.q_eff_fix[irbx, irby, irbz, rbi ] / dV
            end
        end
    end
    return ρ
end
