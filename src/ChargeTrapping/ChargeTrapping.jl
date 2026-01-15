abstract type AbstractChargeTrappingModel{T <: SSDFloat} end

function _calculate_signal( 
        ::AbstractChargeTrappingModel{T},
        path::AbstractVector{CartesianPoint{T}}, 
        pathtimestamps::AbstractVector{T}, 
        charge::T,          
        wpot::Interpolations.Extrapolation{T, 3}, 
        point_types::PointTypes{T, N, S}
    ) where {T <: SSDFloat, N, S}
    throw("For the chosen charge trapping model, no method for `_calculate_signal` is implemented.")
end


"""
    struct NoChargeTrappingModel{T <: SSDFloat} <: AbstractChargeTrappingModel{T}
        
Charge trapping model, in which no charges are trapped during the charge drift.

This model is the default when no charge trapping model is defined in the configuration file.
"""
struct NoChargeTrappingModel{T <: SSDFloat} <: AbstractChargeTrappingModel{T} end

function _signal!(
    ::NoChargeTrappingModel{T},
    q::Base.RefValue{T},
    w::T,
    ::Base.RefValue{T};
    kwargs...
) where {T <: SSDFloat}
    
    return w * q[]
end

function _calculate_signal( 
        ctm::NoChargeTrappingModel{T},
        path::AbstractVector{CartesianPoint{T}}, 
        pathtimestamps::AbstractVector{T}, 
        charge::T,          
        wpot::Interpolations.Extrapolation{T, 3}, 
        point_types::PointTypes{T, N, S}
    )::Vector{T} where {T <: SSDFloat, N, S}

    tmp_signal::Vector{T} = Vector{T}(undef, length(pathtimestamps))

    q = Ref(charge)
    running_sum = Ref(zero(T))
    
    @inbounds for i in eachindex(tmp_signal)
        w = get_interpolation(wpot, path[i], S)::T
        tmp_signal[i] = _signal!(ctm, q, w, running_sum)
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

function _signal!(
    ::BoggsChargeTrappingModel{T},
    q::Base.RefValue{T},
    w::T,
    running_sum::Base.RefValue{T};
    nσ::T=zero(T),
    Δl::T=zero(T),
    kwargs...
) where {T <: SSDFloat}

    Δq::T = q[] * nσ * Δl
    q[] -= Δq
    running_sum[] += w * Δq
    return running_sum[] + w * q[]
end

function _calculate_signal( 
        ctm::BoggsChargeTrappingModel{T},
        path::AbstractVector{CartesianPoint{T}}, 
        pathtimestamps::AbstractVector{T}, 
        charge::T,          
        wpot::Interpolations.Extrapolation{T, 3},
        point_types::PointTypes{T, N, S}
    )::Vector{T} where {T <: SSDFloat, N, S}

    tmp_signal::Vector{T} = Vector{T}(undef, length(pathtimestamps))

    q = Ref(charge)
    running_sum = Ref(zero(T))
    
    vth::T = sqrt(3 * kB * ctm.temperature / (ifelse(charge > 0, ctm.meffh, ctm.meffe) * me)) # in m/s
    nσ::T = ifelse(charge > 0, ctm.nσh, ctm.nσe)
    
    @inbounds for i in eachindex(tmp_signal)
        Δldrift::T = (i > 1) ? norm(path[i] - path[i-1]) : zero(T)
        Δl::T = Δldrift > 0 ? hypot(Δldrift, vth * (pathtimestamps[i] - pathtimestamps[i-1])) : zero(T)
        w::T = i > 1 ? get_interpolation(wpot, path[i], S) : zero(T)

        tmp_signal[i] = _signal!( ctm, q, w, running_sum; nσ=nσ, Δl=Δl )
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

## Fields
* `τh::T`: Lifetime for holes in bulk volume (sensitive region) (default: `τh = 1ms`).
* `τe::T`: Lifetime for electrons in bulk volume (default: `τe = 1ms`).
> Note: All τ must be much bigger than `Δt`.

See also [Charge Trapping Models](@ref).
"""

struct ConstantLifetimeChargeTrappingModel{T <: SSDFloat} <: AbstractChargeTrappingModel{T}
    τh::T
    τe::T
end

function check_lifetime(model::AbstractChargeTrappingModel{T}, charge::T, pathtimestamps::AbstractVector{T}) where {T <: SSDFloat}
    τ = ifelse(charge > 0, model.τh, model.τe)
    Δt_minimum = minimum(diff(pathtimestamps))
    if τ < Δt_minimum
        throw(ArgumentError("The carrier lifetime should be at least bigger than Δt"))
    end
    return τ
end

function _signal!(
    ::ConstantLifetimeChargeTrappingModel{T},
    q::Base.RefValue{T},
    w::T,
    running_sum::Base.RefValue{T};
    τ::T=T(Inf),
    Δt::T=zero(T),
    kwargs...
) where {T <: SSDFloat}

    Δq::T = q[] * Δt / τ
    q[] -= Δq
    running_sum[] += w * Δq
    return running_sum[] + w * q[]
end

function _calculate_signal( 
        ctm::ConstantLifetimeChargeTrappingModel{T},
        path::AbstractVector{CartesianPoint{T}}, 
        pathtimestamps::AbstractVector{T}, 
        charge::T,          
        wpot::Interpolations.Extrapolation{T, 3},
        point_types::PointTypes{T, N, S}
    )::Vector{T} where {T <: SSDFloat, N, S}
    
    tmp_signal::Vector{T} = Vector{T}(undef, length(pathtimestamps))

    q = Ref(charge)
    running_sum = Ref(zero(T))

    τ::T = check_lifetime(ctm, charge, pathtimestamps)
    
    @inbounds for i in eachindex(tmp_signal)
        
        Δt::T = (i > 1) ? (pathtimestamps[i] - pathtimestamps[i-1]) : zero(T)
        w::T = i > 1 ? get_interpolation(wpot, path[i], S) : zero(T)
        
        tmp_signal[i] = _signal!( ctm, q, w, running_sum; τ=τ, Δt=Δt )
    end
 
    tmp_signal
end

ConstantLifetimeChargeTrappingModel(args...; T::Type{<:SSDFloat}, kwargs...) = ConstantLifetimeChargeTrappingModel{T}(args...; kwargs...)
function ConstantLifetimeChargeTrappingModel{T}(config_dict::AbstractDict = Dict()) where {T <: SSDFloat}
    τh::T = ustrip(u"s", 1u"ms")
    τe::T = ustrip(u"s", 1u"ms")

    if haskey(config_dict, "model") && !haskey(config_dict, "parameters")
        throw(ConfigFileError("`ConstantLifetimeChargeTrappingModel` does not have `parameters`"))
    end

    parameters = haskey(config_dict, "parameters") ? config_dict["parameters"] : config_dict

    allowed_keys = ("τh", "τe")
    k = filter(k -> !(k in allowed_keys), keys(parameters))
    !isempty(k) && @warn "The following keys will be ignored: $(k).\nAllowed keys are: $(allowed_keys)"

    if haskey(parameters, "τh")   τh =     _parse_value(T, parameters["τh"], internal_time_unit) end
    if haskey(parameters, "τe")   τe =     _parse_value(T, parameters["τe"], internal_time_unit) end
    ConstantLifetimeChargeTrappingModel{T}(τh, τe)
end


"""
    struct CombinedChargeTrappingModel{T <: SSDFloat} <: AbstractChargeTrappingModel{T}
        
Combined Charge Trapping Model (CTM) for the bulk and inactive layer

## Fields
* `bulk_charge_trapping_model::AbstractChargeTrappingModel{T}`: CTM inside the bulk.
* `inactive_charge_trapping_model::AbstractChargeTrappingModel{T}`: CTM inside the inactive layer.
* `inactive_layer_geometry::G`: Inactive layer geometry.

See also [Charge Trapping Models](@ref).
"""

_ctmodel_name(::BoggsChargeTrappingModel) = "Boggs"
_ctmodel_name(::ConstantLifetimeChargeTrappingModel) = "ConstantLifetime"
_ctmodel_name(::NoChargeTrappingModel) = "NoChargeTrapping"

has_inactive_layer(point_types::Array{PointType, 3}) = any(is_in_inactive_layer, point_types)

struct CombinedChargeTrappingModel{T <: SSDFloat, BCTM <: AbstractChargeTrappingModel{T}, ICTM <: AbstractChargeTrappingModel{T}, G <: Union{<:AbstractGeometry, Nothing}} <: AbstractChargeTrappingModel{T}
    bulk_charge_trapping_model::BCTM
    inactive_charge_trapping_model::ICTM
    inactive_layer_geometry::G
end


function _calculate_signal( 
        ctm::CombinedChargeTrappingModel{T},
        path::AbstractVector{CartesianPoint{T}}, 
        pathtimestamps::AbstractVector{T}, 
        charge::T,          
        wpot::Interpolations.Extrapolation{T, 3},
        point_types::PointTypes{T, N, S}
    )::Vector{T} where {T <: SSDFloat, N, S}
    
    tmp_signal::Vector{T} = Vector{T}(undef, length(pathtimestamps))

    q = Ref(charge)
    running_sum = Ref(zero(T))

    inactive_layer_exist = !isnothing(ctm.inactive_layer_geometry) || has_inactive_layer(point_types.data)

    bulk_ctmodel = _ctmodel_name(ctm.bulk_charge_trapping_model)
    inactive_layer_ctmodel = _ctmodel_name(ctm.inactive_charge_trapping_model)

    τ::T = zero(T)
    τ_inactive::T = zero(T)
    
    vth::T = zero(T)
    nσ::T = zero(T)

    if bulk_ctmodel == "ConstantLifetime"
        τ = check_lifetime(ctm.bulk_charge_trapping_model, charge, pathtimestamps)
    elseif bulk_ctmodel == "Boggs"
        vth = sqrt(3 * kB * ctm.bulk_charge_trapping_model.temperature / (ifelse(charge > 0, ctm.bulk_charge_trapping_model.meffh, ctm.bulk_charge_trapping_model.meffe) * me)) # in m/s
        nσ = ifelse(charge > 0, ctm.bulk_charge_trapping_model.nσh, ctm.bulk_charge_trapping_model.nσe)
    end

    if inactive_layer_ctmodel == "ConstantLifetime"
        τ_inactive = check_lifetime(ctm.inactive_charge_trapping_model, charge, pathtimestamps)
    end

    
    @inbounds for i in eachindex(tmp_signal)

        in_inactive_region = in_inactive_layer(path[i], ctm.inactive_layer_geometry, point_types)
        
        Δt::T = (i > 1) ? (pathtimestamps[i] - pathtimestamps[i-1]) : zero(T)
        w::T = i > 1 ? get_interpolation(wpot, path[i], S) : zero(T)

        Δldrift::T = (i > 1) ? norm(path[i] - path[i-1]) : zero(T)
        Δl::T = Δldrift > 0 ? hypot(Δldrift, vth * (pathtimestamps[i] - pathtimestamps[i-1])) : zero(T)
        
        if in_inactive_region
            tmp_signal[i] = _signal!(ctm.inactive_charge_trapping_model, q, w, running_sum;
                                     τ=τ_inactive, Δt=Δt)
        else
            tmp_signal[i] = _signal!(ctm.bulk_charge_trapping_model, q, w, running_sum;
                                     nσ=nσ, Δl=Δl, τ=τ, Δt=Δt)
        end
    end
    tmp_signal
end



CombinedChargeTrappingModel(args...; T::Type{<:SSDFloat}, kwargs...) = CombinedChargeTrappingModel{T}(args...; kwargs...)
function CombinedChargeTrappingModel{T}(config_dict::AbstractDict = Dict(); temperature::RealQuantity = T(78)) where {T <: SSDFloat}

    if haskey(config_dict, "model") && config_dict["model"] !== nothing && !(haskey(config_dict, "parameters"))
        throw(ConfigFileError("`CombinedChargeTrappingModel` does not have `parameters`"))
    end
    if haskey(config_dict, "model_inactive") && config_dict["model_inactive"] !== nothing && !(haskey(config_dict, "parameters_inactive"))
        throw(ConfigFileError("`CombinedChargeTrappingModel` does not have `parameters_inactive`"))
    end

    bulk_charge_trapping_model = begin
        if !haskey(config_dict, "model")
            NoChargeTrappingModel{T}()
        else
            model = config_dict["model"]
            if isnothing(model)
                throw(ConfigFileError("`model` is defined but empty in config. Remove it or provide a valid name."))
        elseif model == "Boggs"
                BoggsChargeTrappingModel{T}(config_dict; temperature)
            elseif model == "ConstantLifetime"
                ConstantLifetimeChargeTrappingModel{T}(config_dict)
            else
                throw(ConfigFileError("Unknown bulk charge trapping model: $model"))
            end
        end
    end
    
    # Rename 'parameters_inactive' to be 'parameters'
    tmp = copy(config_dict)
    if haskey(tmp, "parameters_inactive")
        if haskey(tmp, "parameters")
            tmp["parameters_bulk"] = tmp["parameters"]
        end
        tmp["parameters"] = tmp["parameters_inactive"]
    end
    
    inactive_charge_trapping_model = if haskey(tmp, "model_inactive") && tmp["model_inactive"] == "ConstantLifetime"     
        ConstantLifetimeChargeTrappingModel{T}(tmp)
    else
        if !isnothing(config_dict["model_inactive"]) @warn "Unknown inactive charge trapping model, running `NoChargeTrappingModel` inside inactive" end
        NoChargeTrappingModel{T}()
    end
    
    inactive_layer_geometry = haskey(config_dict, "inactive_layer_geometry") ? config_dict["inactive_layer_geometry"] : nothing
    
CombinedChargeTrappingModel{T, typeof(bulk_charge_trapping_model), typeof(inactive_charge_trapping_model), typeof(inactive_layer_geometry)}(
    bulk_charge_trapping_model,
    inactive_charge_trapping_model,
    inactive_layer_geometry
)
end
