abstract type AbstractChargeTrappingModel{T <: SSDFloat} end

function _calculate_signal( 
        ::AbstractChargeTrappingModel{T},
        path::AbstractVector{CartesianPoint{T}}, 
        pathtimestamps::AbstractVector{T}, 
        charge::T,          
        wpot::Interpolations.Extrapolation{T, 3}, 
        S::CoordinateSystemType
    ) where {T <: SSDFloat}
    throw("For the chosen charge trapping model, no method for `_calculate_signal` is implemented.")
end


"""
    struct NoChargeTrappingModel{T <: SSDFloat} <: AbstractChargeTrappingModel{T}
        
Charge trapping model, in which no charges are trapped during the charge drift.

This model is the default when no charge trapping model is defined in the configuration file.
"""
struct NoChargeTrappingModel{T <: SSDFloat} <: AbstractChargeTrappingModel{T} end

function _calculate_signal( 
        ::NoChargeTrappingModel{T}, 
        path::AbstractVector{CartesianPoint{T}}, 
        pathtimestamps::AbstractVector{T}, 
        charge::T,          
        wpot::Interpolations.Extrapolation{T, 3}, 
        S::CoordinateSystemType
    )::Vector{T} where {T <: SSDFloat}

    tmp_signal::Vector{T} = Vector{T}(undef, length(pathtimestamps))
    @inbounds for i in eachindex(tmp_signal)
        tmp_signal[i] = get_interpolation(wpot, path[i], S)::T * charge
    end

    tmp_signal
end


"""
    struct BoggsChargeTrappingModel{T <: SSDFloat} <: AbstractChargeTrappingModel{T}
        
Charge trapping model presented in [Steve Boggs (2023)](https://doi.org/10.1016/j.nima.2023.168756).

## Fields
* `nσe::T`: Trapping product for electrons.
* `nσh::T`: Trapping product for holes.
* `Temp::T`: Temperature of the crystal.

"""
Parameters.@with_kw struct BoggsChargeTrappingModel{T <: SSDFloat} <: AbstractChargeTrappingModel{T} 
    nσe::T = ustrip(u"m^-1", inv(1020u"cm"))
    nσh::T = ustrip(u"m^-1", inv(2040u"cm"))
    Temp::T = T(80)
end

function _calculate_signal( 
        ctm::BoggsChargeTrappingModel{T}, 
        path::AbstractVector{CartesianPoint{T}}, 
        pathtimestamps::AbstractVector{T}, 
        charge::T,          
        wpot::Interpolations.Extrapolation{T, 3}, 
        S::CoordinateSystemType
    )::Vector{T} where {T <: SSDFloat}

    vth::T = sqrt(3 * kB * ctm.Temp / (ifelse(charge > 0, 0.21, 0.12) * me)) # in m/s
    nσ::T = ifelse(charge > 0, ctm.nσh, ctm.nσe)
    tmp_signal::Vector{T} = Vector{T}(undef, length(pathtimestamps))

    running_sum::T = zero(T)
    q::T = charge
    @inbounds for i in eachindex(tmp_signal)
        Δldrift::T = (i > 1) ? norm(path[i] .- path[i-1]) : zero(T)
        Δl::T = Δldrift > 0 ? hypot(Δldrift, vth * (pathtimestamps[i] - pathtimestamps[i-1])) : zero(T)
        Δq::T = q * nσ * Δl
        q -= Δq
        w::T = i > 1 ? get_interpolation(wpot, path[i], S) : zero(T)
        running_sum += w * Δq
        tmp_signal[i] = running_sum + w * q
    end

    tmp_signal
end



