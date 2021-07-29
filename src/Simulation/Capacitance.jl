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
    r0_handling::Bool = typeof(ef.grid.axes[1]).parameters[2] == :r0

    ax1::Vector{T} = collect(ef.grid[1])
    ax2::Vector{T} = collect(ef.grid[2])
    ax3::Vector{T} = collect(ef.grid[3])
    mp1::Vector{T} = midpoints(get_extended_ticks(ef.grid[1]))
    mp2::Vector{T} = midpoints(get_extended_ticks(ef.grid[2]))
    mp3::Vector{T} = midpoints(get_extended_ticks(ef.grid[3]))
    Δax1::Vector{T} = diff(ax1)
    Δax2::Vector{T} = diff(ax2)
    Δax3::Vector{T} = diff(ax3)
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