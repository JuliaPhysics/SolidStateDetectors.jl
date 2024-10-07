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

NoChargeTrappingModel(args...; T::Type{<:SSDFloat} = Float32, kwargs...) = NoChargeTrappingModel{T}(args...; kwargs...)
NoChargeTrappingModel{T}(config_dict::AbstractDict; kwargs...) where {T <: SSDFloat} = NoChargeTrappingModel{T}()


"""
    struct BoggsChargeTrappingModel{T <: SSDFloat} <: AbstractChargeTrappingModel{T}
        
Charge trapping model presented in [Boggs _et al._ (2023)](https://doi.org/10.1016/j.nima.2023.168756).

## Fields
* `nσe::T`: Trapping product for electrons (default: `(nσe)^-1 = 1020cm`).
* `nσh::T`: Trapping product for holes (default: `(nσh)^-1 = 2040cm`).
* `temperature::T`: Temperature of the crystal (default: `78K`).

See also [Charge Trapping Models](@ref).
"""
struct BoggsChargeTrappingModel{T <: SSDFloat} <: AbstractChargeTrappingModel{T} 
    nσe::T  # in m^-1
    nσh::T  # in m^-1
    meffe::T # in units of me
    meffh::T # in units of me
    temperature::T # in K
end

function _calculate_signal( 
        ctm::BoggsChargeTrappingModel{T}, 
        path::AbstractVector{CartesianPoint{T}}, 
        pathtimestamps::AbstractVector{T}, 
        charge::T,          
        wpot::Interpolations.Extrapolation{T, 3}, 
        S::CoordinateSystemType
    )::Vector{T} where {T <: SSDFloat}

    vth::T = sqrt(3 * kB * ctm.temperature / (ifelse(charge > 0, ctm.meffh, ctm.meffe) * me)) # in m/s
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


BoggsChargeTrappingModel(args...; T::Type{<:SSDFloat}, kwargs...) = BoggsChargeTrappingModel{T}(args...; kwargs...)
function BoggsChargeTrappingModel{T}(config_dict::AbstractDict = Dict(); temperature::RealQuantity = T(78)) where {T <: SSDFloat}
    nσe::T = ustrip(u"m^-1", inv(1020u"cm"))
    nσh::T = ustrip(u"m^-1", inv(2040u"cm"))
    meffe::T = 0.12
    meffh::T = 0.21
    temperature::T = _parse_value(T, temperature, internal_temperature_unit)
    if haskey(config_dict, "parameters")
        if haskey(config_dict, "nσe")   nσe =     _parse_value(T, config_dict["nσe"], internal_length_unit^-1) end
        if haskey(config_dict, "nσe-1") nσe = inv(_parse_value(T, config_dict["nσe-1"], internal_length_unit)) end
        if haskey(config_dict, "nσh")   nσh =     _parse_value(T, config_dict["nσh"], internal_length_unit^-1) end
        if haskey(config_dict, "nσh-1") nσh = inv(_parse_value(T, config_dict["nσh-1"], internal_length_unit)) end
        if haskey(config_dict, "meffe") meffe =   _parse_value(T, config_dict["meffe"], Unitful.NoUnits) end
        if haskey(config_dict, "meffh") meffh =   _parse_value(T, config_dict["meffh"], Unitful.NoUnits) end
        if haskey(config_dict, "temperature") temperature = _parse_value(T, config_dict["temperature"], internal_temperature_unit) end
    end
    BoggsChargeTrappingModel{T}(nσe, nσh, meffe, meffh, temperature)
end



