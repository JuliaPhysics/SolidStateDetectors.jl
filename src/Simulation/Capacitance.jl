function _get_abs_bias_voltage(det::SolidStateDetector{T}) where {T <: SSDFloat}
    potentials::Vector{T} = map(c -> c.potential, det.contacts)
    return (maximum(potentials) - minimum(potentials)) * u"V"
end

"""
    calculate_stored_energy(sim::Simulation{T}) where {T <: SSDFloat}

Calculates and returns the energy stored in the [`ElectricField`](@ref) of a 
[`SolidStateDetector`](@ref) in a given [`Simulation`](@ref) in units of J.

## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) with `sim.detector` for which the stored energy is calculated.

!!! note 
    This method only works if `sim.electric_field` has already been calculated and is not `missing`.
"""
function calculate_stored_energy(sim::Simulation{T}) where {T <: SSDFloat}
    @assert !ismissing(sim.electric_field) "Electric field has not been calculated yet. Please run `calculate_electric_field!(sim)` first."
    calculate_stored_energy(sim.electric_field, sim.ϵ_r)
end

function calculate_stored_energy(ef::ElectricField{T,3,S}, ϵ::DielectricDistribution{T,3,S}) where {T <: SSDFloat, S}
    W::T = 0

    cylindric::Bool = S == Cylindrical
    cartesian::Bool = !cylindric

    ax1::Vector{T} = collect(ef.grid[1])
    ax2::Vector{T} = collect(ef.grid[2])
    ax3::Vector{T} = collect(ef.grid[3])
    mp1::Vector{T} = midpoints(get_extended_ticks(ef.grid[1]))
    mp2::Vector{T} = midpoints(get_extended_ticks(ef.grid[2]))
    mp3::Vector{T} = midpoints(get_extended_ticks(ef.grid[3]))
    Δmp1::Vector{T} = diff(mp1)
    Δmp2::Vector{T} = diff(mp2)
    Δmp3::Vector{T} = diff(mp3)

    w1r::Vector{T} = inv.(Δmp1) .* (mp1[2:end] .- ax1)
    w1l::Vector{T} = inv.(Δmp1) .* (ax1 - mp1[1:end-1])
    w2r::Vector{T} = inv.(Δmp2) .* (mp2[2:end] .- ax2)
    w2l::Vector{T} = inv.(Δmp2) .* (ax2 - mp2[1:end-1])
    w3r::Vector{T} = inv.(Δmp3) .* (mp3[2:end] .- ax3)
    w3l::Vector{T} = inv.(Δmp3) .* (ax3 - mp3[1:end-1])

    if cylindric
        mp1[1] = 0
        mp1[end] = ax1[end]
        Δmp1 = ((mp1[2:end].^2) .- (mp1[1:end-1].^2)) ./ 2
    end
    V::T = 0
    for i3 in 1:size(ϵ, 3)-1
        _Δmp3::T = Δmp3[i3]
        if (i3 == 1 || i3 == size(ϵ, 3)-1) _Δmp3 /= 2 end
        for i2 in 1:size(ϵ, 2)-1
            _Δmp2::T = Δmp2[i2]
            if (cartesian && (i2 == 1 || i2 == size(ϵ, 2)-1)) _Δmp2 /= 2 end
            for i1 in 1:size(ϵ, 1)-1
                _Δmp1::T = Δmp1[i1]
                if (cartesian && (i1 == 1 || i1 == size(ϵ, 1)-1)) _Δmp1 /= 2 end
                ev::SArray{Tuple{3},Float32,1,3} = ef.data[i1, i2, i3]
                dV::T = _Δmp3 * _Δmp2 * _Δmp1
                _ϵ::T = sum([
                    ϵ[i1, i2, i3]             * w1l[i1] * w2l[i2] * w3l[i3],
                    ϵ[i1 + 1, i2, i3]         * w1r[i1] * w2l[i2] * w3l[i3],
                    ϵ[i1, i2 + 1, i3]         * w1l[i1] * w2r[i2] * w3l[i3],
                    ϵ[i1, i2, i3 + 1]         * w1l[i1] * w2l[i2] * w3r[i3],
                    ϵ[i1 + 1, i2 + 1, i3]     * w1r[i1] * w2r[i2] * w3l[i3],
                    ϵ[i1 + 1, i2, i3 + 1]     * w1r[i1] * w2l[i2] * w3r[i3],
                    ϵ[i1, i2 + 1, i3 + 1]     * w1l[i1] * w2r[i2] * w3r[i3],
                    ϵ[i1 + 1, i2 + 1, i3 + 1] * w1r[i1] * w2r[i2] * w3r[i3]
                ])
                V += dV
                W += sum(ev.^2) * dV * _ϵ
            end
        end
    end
    E =  W * ϵ0 / 2 * u"J"
    if cylindric && (size(ϵ, 2) - 1 != size(ef, 2))
        E *= size(ef, 2) / (size(ϵ, 2) - 1)
    end
    return E
end

export calculate_capacitance
"""
    calculate_capacitance(sim::Simulation{T}) where {T <: SSDFloat}

Calculates and returns the capacitance of a [`SolidStateDetector`](@ref) in a given [`Simulation`](@ref) in units of pF.

## Arguments
* `sim::Simulation{T}`: [`Simulation`](@ref) with `sim.detector` for which the capacitance is calculated.

!!! note 
    This method only works if `sim.electric_field` has already been calculated and is not `missing`.
"""
function calculate_capacitance(sim::Simulation{T}) where {T <: SSDFloat}
    @assert !ismissing(sim.electric_field) "Electric field has not been calculated yet. Please run `calculate_electric_field!(sim)` first."
    W = calculate_stored_energy(sim)
    return uconvert(u"pF", 2 * W / (_get_abs_bias_voltage(sim.detector)^2))
end


function calculate_capacitance(sim::Simulation{T}, ::Type{ElectricPotential}; consider_multiplicity::Bool = true) where {T <: SSDFloat}
    @assert !ismissing(sim.electric_potential) "Electric potential has not been calculated yet. Please run `calculate_electric_potential!(sim)` first."
    W = calculate_stored_energy(sim, ElectricPotential; consider_multiplicity)
    return uconvert(u"pF", 2 * W / (_get_abs_bias_voltage(sim.detector)^2))
end

function calculate_stored_energy(sim::Simulation{T}, ::Type{ElectricPotential}; consider_multiplicity::Bool = true) where {T <: SSDFloat}
    calculate_stored_energy(sim.electric_potential, sim.ϵ_r; consider_multiplicity)
end

function calculate_stored_energy(ep::ScalarPotential{T,3,CS}, ϵ::DielectricDistribution{T,3,CS}; consider_multiplicity::Bool = true) where {T <: SSDFloat, CS}
    cylindrical = CS == Cylindrical
    phi_2D = cylindrical && size(ep, 2) == 1
    ep3d = phi_2D ? get_2π_potential(ep, n_points_in_φ = 2) : _get_closed_potential(ep)
    grid = ep3d.grid
    W::T = 0
    for i3 in 1:size(grid, 3)-1
        for i2 in 1:size(grid, 2)-1
            for i1 in 1:size(grid, 1)-1
                w1, w2, w3 = voxel_widths(grid, i1, i2, i3)
                dV = voxel_volume(grid, i1, i2, i3, w1, w2, w3)

                _ϵ = ϵ.data[i1 + 1, i2 + 1, i3 + 1]

                ep000 = ep3d.data[i1    , i2    , i3    ]
                ep100 = ep3d.data[i1 + 1, i2    , i3    ]
                ep010 = ep3d.data[i1    , i2 + 1, i3    ]
                ep110 = ep3d.data[i1 + 1, i2 + 1, i3    ]
                ep001 = ep3d.data[i1    , i2    , i3 + 1]
                ep101 = ep3d.data[i1 + 1, i2    , i3 + 1]
                ep011 = ep3d.data[i1    , i2 + 1, i3 + 1]
                ep111 = ep3d.data[i1 + 1, i2 + 1, i3 + 1]

                efv1 = ( (ep100 - ep000) + (ep110 - ep010) + (ep101 - ep001) + (ep111 - ep011) ) / (4 * w1)
                efv2 = if cylindrical
                    _w2 = (grid[2].ticks[i2 + 1] - grid[2].ticks[i2])
                    if i1 == 1
                        ((ep110 - ep100)/(_w2*grid[1].ticks[i1+1]) +
                        (ep111 - ep101)/(_w2*grid[1].ticks[i1+1])) / 2
                    else
                        ((ep010 - ep000)/(_w2*grid[1].ticks[i1]) +
                         (ep110 - ep100)/(_w2*grid[1].ticks[i1+1]) +
                         (ep011 - ep001)/(_w2*grid[1].ticks[i1]) +
                         (ep111 - ep101)/(_w2*grid[1].ticks[i1+1])) / 4
                    end
                else
                    ( (ep010 - ep000) + (ep110 - ep100) + (ep011 - ep001) + (ep111 - ep101) ) / (4 * w2)
                end
                efv3 = ( (ep001 - ep000) + (ep101 - ep100) + (ep011 - ep010) + (ep111 - ep110) ) / (4 * w3)
                W += sum((efv1, efv2, efv3).^2) * dV * _ϵ
            end
        end
    end
    E = W * ϵ0 / 2 * u"J"
    phi_2D && (E *= 2)
    return consider_multiplicity ? E * multiplicity(grid) : E
end

function calculate_capacitance(sim::Simulation{T}, ::Type{WeightingPotential}, contact_id::Int = 1; consider_multiplicity::Bool = true) where {T <: SSDFloat}
    @assert !ismissing(sim.weighting_potentials[contact_id]) "Weighting potential of contact $contact_id has not been calculated yet. Please run `calculate_weighting_potential!(sim, $contact_id)` first."
    W = calculate_stored_energy(sim, WeightingPotential, contact_id; consider_multiplicity)
    return uconvert(u"pF", 2 * W / u"V^2")
end
function calculate_stored_energy(sim::Simulation{T, CS}, ::Type{WeightingPotential}, contact_id::Int = 1; consider_multiplicity::Bool = true) where {T <: SSDFloat, CS}
    ϵ_r = DielectricDistribution(DielectricDistributionArray(PotentialSimulationSetupRB(sim.detector, sim.weighting_potentials[contact_id].grid, sim.medium, sim.weighting_potentials[contact_id].data,
                weighting_potential_contact_id = contact_id, use_nthreads = _guess_optimal_number_of_threads_for_SOR(size(sim.weighting_potentials[contact_id].grid), Base.Threads.nthreads(), CS),    
                not_only_paint_contacts = false, paint_contacts = false)), sim.weighting_potentials[contact_id].grid)
    calculate_stored_energy(sim.weighting_potentials[contact_id], ϵ_r; consider_multiplicity)
end
