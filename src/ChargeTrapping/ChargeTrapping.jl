abstract type AbstractChargeTrappingModel{T <: SSDFloat} end

function _calculate_signal( 
        ::AbstractChargeTrappingModel{T},
        path::AbstractVector{CartesianPoint{T}}, 
        pathtimestamps::AbstractVector{T}, 
        charge::T,          
        wpot::Interpolations.Extrapolation{T, 3}, 
        point_types::Union{PointTypes{T, N, S}, Nothing} = nothing
    ) where {T <: SSDFloat, N, S}
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
        point_types::Union{PointTypes{T, N, S}, Nothing} = nothing
    )::Vector{T} where {T <: SSDFloat, N, S}

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
        point_types::Union{PointTypes{T, N, S}, Nothing} = nothing
    )::Vector{T} where {T <: SSDFloat, N, S}

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

    if haskey(config_dict, "model") && !haskey(config_dict, "parameters")
        throw(ConfigFileError("`BoggsChargeTrappingModel` does not have `parameters`"))
    end

    parameters = haskey(config_dict, "parameters") ? config_dict["parameters"] : config_dict
    
    allowed_keys = ("nσe","nσe-1","nσh","nσh-1","meffe","meffh","temperature")
    k = filter(k -> !(k in allowed_keys), keys(parameters))
    !isempty(k) && @warn "The following keys will be ignored: $(k).\nAllowed keys are: $(allowed_keys)"

    if haskey(parameters, "nσe")   nσe =     _parse_value(T, parameters["nσe"], internal_length_unit^-1) end
    if haskey(parameters, "nσe-1") nσe = inv(_parse_value(T, parameters["nσe-1"], internal_length_unit)) end
    if haskey(parameters, "nσh")   nσh =     _parse_value(T, parameters["nσh"], internal_length_unit^-1) end
    if haskey(parameters, "nσh-1") nσh = inv(_parse_value(T, parameters["nσh-1"], internal_length_unit)) end
    if haskey(parameters, "meffe") meffe =   _parse_value(T, parameters["meffe"], Unitful.NoUnits) end
    if haskey(parameters, "meffh") meffh =   _parse_value(T, parameters["meffh"], Unitful.NoUnits) end
    if haskey(parameters, "temperature") temperature = _parse_value(T, parameters["temperature"], internal_temperature_unit) end
    BoggsChargeTrappingModel{T}(nσe, nσh, meffe, meffh, temperature)
end




"""
    struct ConstantLifetimeChargeTrappingModel{T <: SSDFloat} <: AbstractChargeTrappingModel{T}
        
This constant-lifetime-based charge trapping model is similar to the Boggs model, which is constant-mean-free-path based.

The model implemented has two sets of parameters for the sensitive (bulk) and inactive (dead layer/surface layer) volume respectively.

## Fields
* `τh::T`: Lifetime for holes in bulk volume (sensitive region) (default: `τh = 1ms`).
* `τe::T`: Lifetime for electrons in bulk volume (default: `τe = 1ms`).
* `τh_inactive::T`: Lifetime for holes in surface layer (dead layer/inactive layer) (default: `τh_inactive = 1ms`).
* `τe_inactive::T`: Lifetime for electrons in surface layer (default: `τe_inactive = 1ms`).
* `inactive_layer_geometry`: The geometry of the inactive layer. If it's not defined, it uses point_types to get the inactive layer geometry, else it uses the user-defined inactive_layer_geometry.
> Note: All τ must be much bigger than `Δt`.

See also [Charge Trapping Models](@ref).
"""

struct ConstantLifetimeChargeTrappingModel{T <: SSDFloat, G <: Union{<:AbstractGeometry, Nothing}} <: AbstractChargeTrappingModel{T}
    τh::T
    τe::T
    τh_inactive::T
    τe_inactive::T
    inactive_layer_geometry::G
end

has_inactive_layer(point_types::Array{PointType, 3}) = any(is_in_inactive_layer, point_types)

function get_charge_carrier_lifetime(
    pt::CartesianPoint{T},
    g::AbstractGeometry{T},
    point_types::PointTypes{T},
    τ::T,
    τ_inactive::T
    )where {T}
    in(pt, g) ? τ_inactive : τ 
end

function get_charge_carrier_lifetime(
    pt::CartesianPoint{T},
    ::Nothing,
    point_types::PointTypes{T},
    τ::T,
    τ_inactive::T
    )where {T}
    is_in_inactive_layer(point_types[pt]) ? τ_inactive : τ
end

function _calculate_signal( 
        ctm::ConstantLifetimeChargeTrappingModel{T}, 
        path::AbstractVector{CartesianPoint{T}}, 
        pathtimestamps::AbstractVector{T}, 
        charge::T,          
        wpot::Interpolations.Extrapolation{T, 3},
        point_types::Union{PointTypes{T, N, S}, Nothing} = nothing
    )::Vector{T} where {T <: SSDFloat, N, S}
    
    tmp_signal::Vector{T} = Vector{T}(undef, length(pathtimestamps))

    running_sum::T = zero(T)
    q::T = charge

    inactive_layer_exist = !isnothing(ctm.inactive_layer_geometry) || has_inactive_layer(point_types.data)

    τ::T = ifelse(charge > 0, ctm.τh, ctm.τe)
    τ_inactive::T = ifelse(charge > 0, ctm.τh_inactive, ctm.τe_inactive)
    Δt_minimum = minimum(diff(pathtimestamps))
    if τ<Δt_minimum || τ_inactive<Δt_minimum 
        throw(ArgumentError("The carrier lifetime should be at least bigger than Δt in the `ConstantLifetimeChargeTrappingModel`"))
    end
    
    @inbounds for i in eachindex(tmp_signal)
        Δt::T = (i > 1) ? (pathtimestamps[i] - pathtimestamps[i-1]) : zero(T)

        τi::T = inactive_layer_exist ? get_charge_carrier_lifetime(path[i], ctm.inactive_layer_geometry, point_types, τ, τ_inactive) : τ

        Δq::T = q * Δt/τi
        q -= Δq
        w::T = i > 1 ? get_interpolation(wpot, path[i], S) : zero(T)
        running_sum += w * Δq
        tmp_signal[i] = running_sum + w * q
    end
 
    tmp_signal
end

ConstantLifetimeChargeTrappingModel(args...; T::Type{<:SSDFloat}, kwargs...) = ConstantLifetimeChargeTrappingModel{T}(args...; kwargs...)
function ConstantLifetimeChargeTrappingModel{T}(config_dict::AbstractDict = Dict()) where {T <: SSDFloat}
    τh::T = ustrip(u"s", 1u"ms")
    τe::T = ustrip(u"s", 1u"ms")
    τh_inactive::T = ustrip(u"s", 80u"ns")
    τe_inactive::T = ustrip(u"s", 80u"ns")

    if haskey(config_dict, "model") && !haskey(config_dict, "parameters")
        throw(ConfigFileError("`ConstantLifetimeChargeTrappingModel` does not have `parameters`"))
    end

    parameters = haskey(config_dict, "parameters") ? config_dict["parameters"] : config_dict
    inactive_layer_geometry = haskey(config_dict, "inactive_layer_geometry") ? config_dict["inactive_layer_geometry"] : nothing

    allowed_keys = ("τh", "τe", "τh_inactive", "τe_inactive")
    k = filter(k -> !(k in allowed_keys), keys(parameters))
    !isempty(k) && @warn "The following keys will be ignored: $(k).\nAllowed keys are: $(allowed_keys)"

    if haskey(parameters, "τh")   τh =     _parse_value(T, parameters["τh"], internal_time_unit) end
    if haskey(parameters, "τe")   τe =     _parse_value(T, parameters["τe"], internal_time_unit) end
    if haskey(parameters, "τh_inactive")   τh_inactive =     _parse_value(T, parameters["τh_inactive"], internal_time_unit) end
    if haskey(parameters, "τe_inactive")   τe_inactive =     _parse_value(T, parameters["τe_inactive"], internal_time_unit) end
    ConstantLifetimeChargeTrappingModel{T, typeof(inactive_layer_geometry)}(τh, τe, τh_inactive, τe_inactive, inactive_layer_geometry)
end
