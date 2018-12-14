struct PrecalculatedWeightsCylindricalRedBlack{T<:AbstractFloat} <: PrecalculatedWeights
    point_types::Array{PointType, 4}
    volume_weights::Array{T, 4}
    charge::Array{T, 4}
    epsilon_r::Array{T, 3}
    wr::Array{T, 2}
    wφ::Array{T, 2}
    wz::Array{T, 2}
    r_range::UnitRange{Int64}
    φ_range::UnitRange{Int64}
    z_range::UnitRange{Int64}
    sor_const::Array{T, 1}
    minimum_voltage::T

    PrecalculatedWeightsCylindricalRedBlack{T}(
            point_types::Array{PointType, 4},
            volume_weights::Array{T, 4},
            charge::Array{T, 4},
            epsilon_r::Array{T, 3},
            wr::Array{T, 2},
            wφ::Array{T, 2},
            wz::Array{T, 2},
            r_range::UnitRange{Int64},
            φ_range::UnitRange{Int64},
            z_range::UnitRange{Int64},
            sor_const::Array{T, 1},
            minimum_voltage::T,
            ) where {T<:AbstractFloat} = new{T}(
                point_types,
                volume_weights,
                charge,
                epsilon_r,
                wr,
                wφ,
                wz,
                r_range,
                φ_range,
                z_range,
                sor_const,
                minimum_voltage,
    )
end

function midpoints(a::AbstractArray)::AbstractArray
    @inbounds r = a[1:end-1]
    @simd for i in eachindex(r)
        @inbounds r[i] += 0.5 * (a[i + 1] - a[i])
    end
    return r
end

function PrecalculatedWeightsCylindricalRedBlack(detector::SolidStateDetector, grid::CylindricalGrid; sor_consts::Vector{<:Real}=[1.4, 1.85], only_2d::Bool=false)::PrecalculatedWeightsCylindricalRedBlack
    T = get_precision_type(detector)

    array_size = size(grid.potential, 1), size(grid.potential, 2), div(size(grid.potential, 3), 2) + mod(size(grid.potential, 3), 2)

    point_types = zeros(PointType, array_size..., 2)
    volume_weights = Array{T, 4}(undef, array_size..., 2)
    charge = Array{T, 4}(undef, array_size..., 2)

    pos_r::T = 0
    pos_φ::T = 0
    pos_z::T = 0

    detector_material_ϵ_r::T = detector.material_detector.ϵ_r
    environment_material_ϵ_r::T = detector.material_environment.ϵ_r

    bias_voltage::T = detector.segment_bias_voltages[1]
    sor_slope = (sor_consts[2] .- sor_consts[1]) / (size(grid.potential, 1) - 1 )
    sor_const::Array{T, 1} = T[ sor_consts[1] + (i - 1) * sor_slope for i in 1:size(grid.potential, 1)]

    odd::Int  = 1
    even::Int = 2

    @inbounds begin
        nr::Int = length(grid.r)
        nφ::Int = length(grid.φ)
        nz::Int = length(grid.z)

        Δr::Array{T, 1} = diff(grid.r)
        Δφ::Array{T, 1} = diff(grid.φ)
        Δz::Array{T, 1} = diff(grid.z)

        r_inv::Array{T, 1} = inv.(grid.r)

        if grid.r[1] == 0 r_inv[1] = inv(grid.r[2] * 0.5) end

        r_ext = Array{T, 1}(undef, length(grid.r) + 2)
        φ_ext = Array{T, 1}(undef, length(grid.φ) + 2)
        z_ext = Array{T, 1}(undef, length(grid.z) + 2)
        r_ext[2:end-1] = grid.r[:]
        r_ext[1] = grid.r[1] - Δr[1]
        r_ext[end] = grid.r[end] + (grid.r[end] - grid.r[end-1])

        φ_ext[2:end-1] = grid.φ[:]
        φ_ext[1] = !only_2d ? grid.φ[1] - (grid.cyclic - grid.φ[end]) : -2π #-π/2
        φ_ext[end] = !only_2d ? grid.cyclic : 2π#π/2
        z_ext[2:end-1] = grid.z[:]
        z_ext[1] = grid.z[1] - Δz[1]
        z_ext[end] = grid.z[end] + Δz[end]

        Δr_ext::Array{T, 1} = diff(r_ext)
        Δφ_ext::Array{T, 1} = diff(φ_ext)
        Δz_ext::Array{T, 1} = diff(z_ext)

        Δr_ext_inv::Array{T, 1} = inv.(Δr_ext)
        Δφ_ext_inv::Array{T, 1} = inv.(Δφ_ext)
        Δz_ext_inv::Array{T, 1} = inv.(Δz_ext)

        mpr::Array{T, 1} = midpoints(r_ext)
        mpφ::Array{T, 1} = midpoints(φ_ext)
        mpz::Array{T, 1} = midpoints(z_ext)

        mpr_inv::Array{T, 1} = inv.(mpr)
        mpφ_inv::Array{T, 1} = inv.(mpφ)
        mpz_inv::Array{T, 1} = inv.(mpz)

        # differences between mitpoints. These are needed for the volumen element of a Point: dV
        Δmpr = diff(mpr)
        Δmpr[1] *= 0.5
        Δmpφ = diff(mpφ)
        Δmpz = diff(mpz)
        Δmpr_inv = inv.(Δmpr)
        Δmpφ_inv = inv.(Δmpφ)
        Δmpz_inv = inv.(Δmpz)

        Δmpr_squared::Array{T, 1} = T(0.5) .* ((mpr[2:end].^2) .- (mpr[1:end-1].^2))
        if grid.r[1] == 0
            Δmpr_squared[1] = T(0.5) * (mpr[2]^2)
            # mpr[1] = 0 # why?
            Δr_ext_inv[1] = Δr_ext_inv[2]  # 0
        end

        # distances between midpoints and real grid points (half distances -> h), needed for weights wrr, wrl, ... & dV (volume element)
        Δhmprr::Array{T, 1} = mpr[2:end] - grid.r
        Δhmprl::Array{T, 1} = grid.r - mpr[1:end - 1]
        Δhmpφr::Array{T, 1} = mpφ[2:end] - grid.φ
        Δhmpφl::Array{T, 1} = grid.φ - mpφ[1:end - 1]
        Δhmpzr::Array{T, 1} = mpz[2:end] - grid.z
        Δhmpzl::Array{T, 1} = grid.z - mpz[1:end - 1]


        epsilon_r  = Array{T, 3}(undef, length(mpr), length(mpφ), length(mpz))
        charge_tmp = Array{T, 3}(undef, length(mpr), length(mpφ), length(mpz))
        @inbounds for iz in 1:size(epsilon_r, 3)
            pos_z = mpz[iz]
            for iφ in 1:size(epsilon_r, 2)
                pos_φ = mpφ[iφ]

                for ir in 2:size(epsilon_r, 1)
                    pos_r = mpr[ir]

                    if contains(detector, pos_r, pos_φ, pos_z)
                        epsilon_r[ir, iφ, iz]::T = detector_material_ϵ_r
                        charge_tmp[ir, iφ, iz]::T = get_charge_density(detector, pos_r, pos_φ, pos_z) * effective_electron_charge_vacuum
                    else
                        epsilon_r[ir, iφ, iz]::T = environment_material_ϵ_r
                        charge_tmp[ir, iφ, iz]::T = 0
                    end
                end
                ir = 1
                pos_r = mpr[ir]
                if grid.r[1] == 0
                    pos_r = grid.r[2] * 0.5
                end
                if contains(detector, pos_r, pos_φ, pos_z)
                    epsilon_r[ir, iφ, iz]::T = detector_material_ϵ_r
                    charge_tmp[ir, iφ, iz]::T = get_charge_density(detector, pos_r, pos_φ, pos_z) * effective_electron_charge_vacuum
                else
                    epsilon_r[ir, iφ, iz]::T = environment_material_ϵ_r
                    charge_tmp[ir, iφ, iz]::T = 0
                end
            end
        end

        # weights for epislon_r adding
        wrr::Array{T, 1} = Δmpr_inv .* Δhmprr
        wrl::Array{T, 1} = Δmpr_inv .* Δhmprl
        if grid.r[1] == 0
            wrl[1]::T = 1 - wrr[1]
        end
        wφr::Array{T, 1} = Δmpφ_inv .* Δhmpφr
        wφl::Array{T, 1} = Δmpφ_inv .* Δhmpφl
        wzr::Array{T, 1} = Δmpz_inv .* Δhmpzr
        wzl::Array{T, 1} = Δmpz_inv .* Δhmpzl

        # weighted epsilon_r
        wrr_eps::T = 1 # for V_r+
        wrl_eps::T = 1 # for V_r-
        wφr_eps::T = 1 # for V_φ+
        wφl_eps::T = 1 # for V_φ-
        wzr_eps::T = 1 # for V_z+
        wzl_eps::T = 1 # for V_z-

        volume_weight::T = 0
        dV::T = 0
        ρ::T = 0

        r_range = range(2, stop = length(r_ext) - 1)
        φ_range = range(2, stop = length(φ_ext) - 1)
        z_range = range(2, stop = length(z_ext) - 1)
        @inbounds for iz in z_range
            inz = iz - 1
            pos_z = z_ext[iz]
            for iφ in φ_range
                inφ = iφ - 1
                pos_φ = φ_ext[iφ]
                for ir in r_range
                    inr = ir - 1;
                    pos_r = r_ext[ir]
                    # if inr == 1 && inz == 1 && inφ == 1
                    # end

                    evenodd = iseven(inr + inφ + inz) ? even : odd
                    rbinds = get_rb_inds(inr, inφ, inz)
                    rbinds = rbinds[1], rbinds[2], rbinds[3], evenodd
                    pwinds = rbinds[1] - 1, rbinds[2] - 1, rbinds[3] - 1, evenodd

                    ρ = 0

                    if inr > 1
                        ρ += charge_tmp[ ir,  iφ,  iz] * wzr[inz] * wrr[inr] * wφr[inφ]
                        ρ += charge_tmp[ ir,  iφ, inz] * wzl[inz] * wrr[inr] * wφr[inφ]
                        ρ += charge_tmp[ ir, inφ,  iz] * wzr[inz] * wrr[inr] * wφl[inφ]
                        ρ += charge_tmp[ ir, inφ, inz] * wzl[inz] * wrr[inr] * wφl[inφ]

                        ρ += charge_tmp[inr,  iφ,  iz] * wzr[inz] * wrl[inr] * wφr[inφ]
                        ρ += charge_tmp[inr,  iφ, inz] * wzl[inz] * wrl[inr] * wφr[inφ]
                        ρ += charge_tmp[inr, inφ,  iz] * wzr[inz] * wrl[inr] * wφl[inφ]
                        ρ += charge_tmp[inr, inφ, inz] * wzl[inz] * wrl[inr] * wφl[inφ]
                    else
                        ρ += charge_tmp[ ir,  iφ,  iz] * wzr[inz] * wφr[inφ]
                        ρ += charge_tmp[ ir,  iφ, inz] * wzl[inz] * wφr[inφ]
                        ρ += charge_tmp[ ir, inφ,  iz] * wzr[inz] * wφl[inφ]
                        ρ += charge_tmp[ ir, inφ, inz] * wzl[inz] * wφl[inφ]
                    end

                    wrr_eps  = epsilon_r[  ir,  iφ, inz + 1] * wφr[inφ] * wzr[inz]
                    wrr_eps += epsilon_r[  ir, inφ, inz + 1] * wφl[inφ] * wzr[inz]
                    wrr_eps += epsilon_r[  ir,  iφ, inz ]    * wφr[inφ] * wzl[inz]
                    wrr_eps += epsilon_r[  ir, inφ, inz ]    * wφl[inφ] * wzl[inz]
                    # # left weight in r: wrr
                    wrl_eps  = epsilon_r[ inr,  iφ, inz + 1] * wφr[inφ] * wzr[inz]
                    wrl_eps += epsilon_r[ inr, inφ, inz + 1] * wφl[inφ] * wzr[inz]
                    wrl_eps += epsilon_r[ inr,  iφ, inz ]    * wφr[inφ] * wzl[inz]
                    wrl_eps += epsilon_r[ inr, inφ, inz ]    * wφl[inφ] * wzl[inz]
                    # right weight in φ: wφr
                    wφr_eps  = epsilon_r[ inr,  iφ, inz + 1] * wrl[inr] * wzr[inz]
                    wφr_eps += epsilon_r[  ir,  iφ, inz + 1] * wrr[inr] * wzr[inz]
                    wφr_eps += epsilon_r[ inr,  iφ, inz ]    * wrl[inr] * wzl[inz]
                    wφr_eps += epsilon_r[  ir,  iφ, inz ]    * wrr[inr] * wzl[inz]
                    # left weight in φ: wφl
                    wφl_eps  = epsilon_r[ inr, inφ, inz + 1] * wrl[inr] * wzr[inz]
                    wφl_eps += epsilon_r[  ir, inφ, inz + 1] * wrr[inr] * wzr[inz]
                    wφl_eps += epsilon_r[ inr, inφ, inz ]    * wrl[inr] * wzl[inz]
                    wφl_eps += epsilon_r[  ir, inφ, inz ]    * wrr[inr] * wzl[inz]
                    # right weight in z: wzr
                    wzr_eps  = epsilon_r[  ir,  iφ, inz + 1] * wrr[inr] * wφr[inφ]
                    wzr_eps += epsilon_r[  ir, inφ, inz + 1] * wrr[inr] * wφl[inφ]
                    wzr_eps += epsilon_r[ inr,  iφ, inz + 1] * wrl[inr] * wφr[inφ]
                    wzr_eps += epsilon_r[ inr, inφ, inz + 1] * wrl[inr] * wφl[inφ]
                    # left weight in z: wzr
                    wzl_eps  = epsilon_r[  ir,  iφ,    inz ] * wrr[inr] * wφr[inφ]
                    wzl_eps += epsilon_r[  ir, inφ,    inz ] * wrr[inr] * wφl[inφ]
                    wzl_eps += epsilon_r[ inr,  iφ,    inz ] * wrl[inr] * wφr[inφ]
                    wzl_eps += epsilon_r[ inr, inφ,    inz ] * wrl[inr] * wφl[inφ]

                    volume_weight = wrr_eps * Δr_ext_inv[ ir] * mpr[ ir] * Δmpφ[inφ] * Δmpz[inz]
                    if inr != 1
                        volume_weight += wrl_eps * Δr_ext_inv[inr] * mpr[inr] * Δmpφ[inφ] * Δmpz[inz]
                    end
                    volume_weight += wφr_eps * r_inv[inr] * Δφ_ext_inv[ iφ] * Δmpr[inr] * Δmpz[inz]
                    volume_weight += wφl_eps * r_inv[inr] * Δφ_ext_inv[inφ] * Δmpr[inr] * Δmpz[inz]
                    volume_weight += wzr_eps * Δz_ext_inv[ iz] * Δmpφ[inφ] * Δmpr_squared[inr]
                    volume_weight += wzl_eps * Δz_ext_inv[inz] * Δmpφ[inφ] * Δmpr_squared[inr]

                    volume_weights[pwinds...] = inv(volume_weight)

                    dV = Δmpz[inz] * Δmpφ[inφ] * Δmpr_squared[inr]

                    if !is_boundary_point(detector, pos_r, pos_φ, pos_z, grid.r, grid.φ, grid.z)
                        point_types[pwinds...] += update_bit
                    else
                        grid.potential[ inr, inφ, inz ] = get_boundary_value( detector, pos_r, pos_φ, pos_z, grid.r)
                    end

                    if inr == 1 pos_r = r_ext[ir + 1] * 0.5 end
                    inside_detector::Bool = contains(detector, pos_r, pos_φ, pos_z)
                    if inside_detector point_types[pwinds...] += pn_junction_bit end

                    charge[pwinds...] = dV * ρ
                end
            end
        end
    end # @inbounds
    r_range = 2:nr+1
    φ_range = 2:nφ+1
    z_range = 2:div(nz + 1,2)+1

    vmins::Array{T, 1}      = zeros(T, length(r_range))
    tmp_array::Array{T, 1}  = zeros(T, length(r_range))

    wr::Array{T, 2} = zeros(T, 7, length(wrr) + 1)
    wr[1, 1:length(wrr)] = wrr
    wr[2, 1:length(wrr)] = wrl
    wr[3, 1:length(wrr)] = Δmpr
    wr[4, :] = mpr
    wr[5, :] = Δr_ext_inv
    wr[6, 1:length(wrr)] = r_inv
    wr[7, 1:length(wrr)] = Δmpr_squared

    wφ::Array{T, 2} = zeros(T, 4, length(wφr) + 1)
    wφ[1, 1:length(wφr)] = wφr
    wφ[2, 1:length(wφr)] = wφl
    wφ[3, 1:length(wφr)] = Δmpφ
    wφ[4, :] = Δφ_ext_inv

    wz::Array{T, 2} = zeros(T, 4, length(wzr) + 1)
    wz[1, 1:length(wzr)] = wzr
    wz[2, 1:length(wzr)] = wzl
    wz[3, 1:length(wzr)] = Δmpz
    wz[4, :] = Δz_ext_inv

    return PrecalculatedWeightsCylindricalRedBlack{T}(
        point_types,
        volume_weights,
        charge,
        epsilon_r,
        wr,
        wφ,
        wz,
        r_range,
        φ_range,
        z_range,
        sor_const,
        bias_voltage
    )
end

function get_PointTypes(grid::CylindricalGrid, pw::PrecalculatedWeightsCylindricalRedBlack)::Array{PointType, 3}
    point_types = zeros(PointType, size(grid.potential)...)
    odd  = 1
    even = 2
    for iz in eachindex(grid.z)
        for iφ in eachindex(grid.φ)
            for ir in eachindex(grid.r)
                evenodd = iseven(ir + iφ + iz) ? even : odd
                rbinds = get_rb_inds(ir, iφ, iz)
                rbinds = rbinds[1] - 1, rbinds[2] - 1, rbinds[3] - 1, evenodd

                point_types[ir, iφ, iz] = pw.point_types[rbinds...]
            end
        end
    end
    return point_types
end

"""
    get_active_volume(grid::CylindricalGrid, pts::PointTypes{T}) where {T}

Returns an approximation of the active volume of the detector by summing up the cell volumes of
all depleted cells.
"""
function get_active_volume(grid::CylindricalGrid, pts::PointTypes{T}) where {T}
    active_volume::T = 0

    only_2d::Bool = length(pts.φ) == 1

    r_ext::Array{T, 1} = Array{T}(undef, length(pts.r) + 2)
    φ_ext::Array{T, 1} = Array{T}(undef, length(pts.φ) + 2)
    z_ext::Array{T, 1} = Array{T}(undef, length(pts.z) + 2)

    r_ext[2:end-1] = pts.r[:]
    r_ext[1] = pts.r[1] - (pts.r[2] - pts.r[1]) # Δr[1]
    r_ext[end] = pts.r[end] + (pts.r[end] - pts.r[end-1])

    φ_ext[2:end-1] = pts.φ[:]
    φ_ext[1] = !only_2d ? pts.φ[1] - (pts.φ[2] - pts.φ[1]) : -π #-π/2
    φ_ext[end] = !only_2d ? grid.cyclic : π
    z_ext[2:end-1] = pts.z[:]
    z_ext[1] = pts.z[1] - (pts.z[2] - pts.z[1])
    z_ext[end] = pts.z[end] + (pts.z[end] - pts.z[end-1])

    Δr_ext::Array{T, 1} = diff(r_ext)
    Δφ_ext::Array{T, 1} = diff(φ_ext)
    Δz_ext::Array{T, 1} = diff(z_ext)
    mpr::Array{T, 1} = midpoints(r_ext)
    mpφ::Array{T, 1} = midpoints(φ_ext)
    mpz::Array{T, 1} = midpoints(z_ext)
    Δmpφ::Array{T, 1} = diff(mpφ)
    Δmpz::Array{T, 1} = diff(mpz)
    Δmpr_squared::Array{T, 1} = T(0.5) .* ((mpr[2:end].^2) .- (mpr[1:end-1].^2))
    if pts.r[1] == 0
        Δmpr_squared[1] = T(0.5) * (mpr[2]^2)
    end

    for iz in eachindex(pts.z)
        for iφ in eachindex(pts.φ)
            for ir in eachindex(pts.r)
                pt::PointType = pts[ir, iφ, iz]
                if (pt & pn_junction_bit > 0) && (pt & undepleted_bit == 0) && (pt & update_bit > 0)
                    active_volume += Δmpz[iz] * Δmpφ[iφ] * Δmpr_squared[ir]
                end
            end
        end
    end
    if grid.cyclic > 0
        active_volume *= 2π / grid.cyclic
    end

    f::T = 10^6
    return active_volume * f * Unitful.cm * Unitful.cm * Unitful.cm
end

function is_fully_depleted(pts::PointTypes{T})::Bool where {T}
    return !(undepleted_bit in (pts.pointtypes .& undepleted_bit))
end


function get_epsilon_r_and_charge(detector::SolidStateDetector, grid::CylindricalGrid; sor_consts::Vector{<:Real}=[1.4, 1.85], only_2d::Bool=false)
    T = get_precision_type(detector)

    array_size = size(grid.potential, 1), size(grid.potential, 2), div(size(grid.potential, 3), 2) + mod(size(grid.potential, 3), 2)

    point_types = zeros(PointType, array_size..., 2)
    volume_weights = Array{T, 4}(undef, array_size..., 2)
    charge = Array{T, 4}(undef, array_size..., 2)

    pos_r::T = 0
    pos_φ::T = 0
    pos_z::T = 0

    detector_material_ϵ_r::T = detector.material_detector.ϵ_r
    environment_material_ϵ_r::T = detector.material_environment.ϵ_r

    bias_voltage::T = detector.segment_bias_voltages[1]
    sor_slope = (sor_consts[2] .- sor_consts[1]) / (size(grid.potential, 1) - 1 )
    sor_const::Array{T, 1} = T[ sor_consts[1] + (i - 1) * sor_slope for i in 1:size(grid.potential, 1)]

    odd::Int  = 1
    even::Int = 2

    @inbounds begin
        nr::Int = length(grid.r)
        nφ::Int = length(grid.φ)
        nz::Int = length(grid.z)

        Δr::Array{T, 1} = diff(grid.r)
        Δφ::Array{T, 1} = diff(grid.φ)
        Δz::Array{T, 1} = diff(grid.z)

        r_inv::Array{T, 1} = inv.(grid.r)

        if grid.r[1] == 0 r_inv[1] = inv(grid.r[2] * 0.5) end

        r_ext = Array{T, 1}(undef, length(grid.r) + 2)
        φ_ext = Array{T, 1}(undef, length(grid.φ) + 2)
        z_ext = Array{T, 1}(undef, length(grid.z) + 2)
        r_ext[2:end-1] = grid.r[:]
        r_ext[1] = grid.r[1] - Δr[1]
        r_ext[end] = grid.r[end] + (grid.r[end] - grid.r[end-1])

        φ_ext[2:end-1] = grid.φ[:]
        φ_ext[1] = !only_2d ? grid.φ[1] - (grid.cyclic - grid.φ[end]) : -2π #-π/2
        φ_ext[end] = !only_2d ? grid.cyclic : 2π#π/2
        z_ext[2:end-1] = grid.z[:]
        z_ext[1] = grid.z[1] - Δz[1]
        z_ext[end] = grid.z[end] + Δz[end]

        Δr_ext::Array{T, 1} = diff(r_ext)
        Δφ_ext::Array{T, 1} = diff(φ_ext)
        Δz_ext::Array{T, 1} = diff(z_ext)

        Δr_ext_inv::Array{T, 1} = inv.(Δr_ext)
        Δφ_ext_inv::Array{T, 1} = inv.(Δφ_ext)
        Δz_ext_inv::Array{T, 1} = inv.(Δz_ext)

        mpr::Array{T, 1} = midpoints(r_ext)
        mpφ::Array{T, 1} = midpoints(φ_ext)
        mpz::Array{T, 1} = midpoints(z_ext)

        mpr_inv::Array{T, 1} = inv.(mpr)
        mpφ_inv::Array{T, 1} = inv.(mpφ)
        mpz_inv::Array{T, 1} = inv.(mpz)

        # differences between mitpoints. These are needed for the volumen element of a Point: dV
        Δmpr = diff(mpr)
        Δmpr[1] *= 0.5
        Δmpφ = diff(mpφ)
        Δmpz = diff(mpz)
        Δmpr_inv = inv.(Δmpr)
        Δmpφ_inv = inv.(Δmpφ)
        Δmpz_inv = inv.(Δmpz)

        Δmpr_squared::Array{T, 1} = T(0.5) .* ((mpr[2:end].^2) .- (mpr[1:end-1].^2))
        if grid.r[1] == 0
            Δmpr_squared[1] = T(0.5) * (mpr[2]^2)
            # mpr[1] = 0 # why?
            Δr_ext_inv[1] = Δr_ext_inv[2]  # 0
        end

        # distances between midpoints and real grid points (half distances -> h), needed for weights wrr, wrl, ... & dV (volume element)
        Δhmprr::Array{T, 1} = mpr[2:end] - grid.r
        Δhmprl::Array{T, 1} = grid.r - mpr[1:end - 1]
        Δhmpφr::Array{T, 1} = mpφ[2:end] - grid.φ
        Δhmpφl::Array{T, 1} = grid.φ - mpφ[1:end - 1]
        Δhmpzr::Array{T, 1} = mpz[2:end] - grid.z
        Δhmpzl::Array{T, 1} = grid.z - mpz[1:end - 1]


        epsilon_r  = Array{T, 3}(undef, length(mpr), length(mpφ), length(mpz))
        charge_tmp = Array{T, 3}(undef, length(mpr), length(mpφ), length(mpz))
        @inbounds for iz in 1:size(epsilon_r, 3)
            pos_z = mpz[iz]
            for iφ in 1:size(epsilon_r, 2)
                pos_φ = mpφ[iφ]
                while pos_φ < 0 
                    pos_φ += grid.cyclic
                end
                for ir in 2:size(epsilon_r, 1)
                    pos_r = mpr[ir]

                    if contains(detector, pos_r, pos_φ, pos_z)
                        epsilon_r[ir, iφ, iz]::T = detector_material_ϵ_r
                        charge_tmp[ir, iφ, iz]::T = get_charge_density(detector, pos_r, pos_φ, pos_z) * effective_electron_charge_vacuum
                    else
                        epsilon_r[ir, iφ, iz]::T = environment_material_ϵ_r
                        charge_tmp[ir, iφ, iz]::T = 0
                    end
                end
                ir = 1
                pos_r = mpr[ir]
                if grid.r[1] == 0
                    pos_r = grid.r[2] * 0.5
                end
                if contains(detector, pos_r, pos_φ, pos_z)
                    epsilon_r[ir, iφ, iz]::T = detector_material_ϵ_r
                    charge_tmp[ir, iφ, iz]::T = get_charge_density(detector, pos_r, pos_φ, pos_z) * effective_electron_charge_vacuum
                else
                    epsilon_r[ir, iφ, iz]::T = environment_material_ϵ_r
                    charge_tmp[ir, iφ, iz]::T = 0
                end
            end
        end
    end
    return epsilon_r, charge_tmp
end