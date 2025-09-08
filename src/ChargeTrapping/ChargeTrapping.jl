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

function _nochargetrapping_signal!(
    w::T,
    charge::T
) where {T <: SSDFloat}

    return w * charge
end

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
        w = get_interpolation(wpot, path[i], S)::T
        tmp_signal[i] = _nochargetrapping_signal!(w, charge)
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

function _Boggs_signal!(
    running_sum::Base.RefValue{T},
    nσ::T,
    Δl::T,
    q::Base.RefValue{T},
    w::T,
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
        point_types::Union{PointTypes{T, N, S}, Nothing} = nothing
    )::Vector{T} where {T <: SSDFloat, N, S}

    tmp_signal::Vector{T} = Vector{T}(undef, length(pathtimestamps))

    q = Ref(charge)
    running_sum = Ref(zero(T))
    
    vth::T = sqrt(3 * kB * ctm.temperature / (ifelse(charge > 0, ctm.meffh, ctm.meffe) * me)) # in m/s
    nσ::T = ifelse(charge > 0, ctm.nσh, ctm.nσe)
    
    @inbounds for i in eachindex(tmp_signal)

        Δldrift::T = (i > 1) ? norm(path[i] .- path[i-1]) : zero(T)
        Δl::T = Δldrift > 0 ? hypot(Δldrift, vth * (pathtimestamps[i] - pathtimestamps[i-1])) : zero(T)
        w::T = i > 1 ? get_interpolation(wpot, path[i], S) : zero(T)

        tmp_signal[i] = _Boggs_signal!( running_sum, nσ, Δl, q, w )
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

function _constantlifetime_signal!(
    running_sum::Base.RefValue{T},
    τ::T,
    Δt::T,
    q::Base.RefValue{T},
    w::T,
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
        point_types::Union{PointTypes{T, N, S}, Nothing} = nothing
    )::Vector{T} where {T <: SSDFloat, N, S}
    
    tmp_signal::Vector{T} = Vector{T}(undef, length(pathtimestamps))

    q = Ref(charge)
    running_sum = Ref(zero(T))

    τ::T = ifelse(charge > 0, ctm.τh, ctm.τe)
    Δt_minimum = minimum(diff(pathtimestamps))
    if τ<Δt_minimum
        throw(ArgumentError("The carrier lifetime should be at least bigger than Δt"))
    end
   
    
    @inbounds for i in eachindex(tmp_signal)
        
        Δt::T = (i > 1) ? (pathtimestamps[i] - pathtimestamps[i-1]) : zero(T)
        w::T = i > 1 ? get_interpolation(wpot, path[i], S) : zero(T)
        
        tmp_signal[i] = _constantlifetime_signal!( running_sum, τ, Δt, q, w )
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

Combinations supported: 
Bulk: Boggs CTM, inactive layer: Constant lifetime CTM
Bulk: Boggs CTM, inactive layer: No CTM
Bulk: No CTM, inactive layer: Constant lifetime CTM
Bulk: No CTM, inactive layer: No CTM
Bulk: Constant lifetime CTM, inactive layer: Constant lifetime CTM
Bulk: Constant lifetime CTM, inactive layer: No CTM

See also [Charge Trapping Models](@ref).
"""

_ctmodel_name(::BoggsChargeTrappingModel) = "Boggs"
_ctmodel_name(::ConstantLifetimeChargeTrappingModel) = "ConstantLifetime"
_ctmodel_name(::NoChargeTrappingModel) = "NoChargeTrapping"
_ctmodel_name(::AbstractChargeTrappingModel) = "Unknown"


has_inactive_layer(point_types::Array{PointType, 3}) = any(is_in_inactive_layer, point_types)

function in_inactive_layer(
    pt::CartesianPoint{T},
    g::AbstractGeometry{T},
    point_types::PointTypes{T}
    )where {T}
    in(pt, g)
end

function in_inactive_layer(
    pt::CartesianPoint{T},
    ::Nothing,
    point_types::PointTypes{T}
    )where {T}
    is_in_inactive_layer(point_types[pt])
end

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
        point_types::Union{PointTypes{T, N, S}, Nothing} = nothing
    )::Vector{T} where {T <: SSDFloat, N, S}
    
    tmp_signal::Vector{T} = Vector{T}(undef, length(pathtimestamps))

    q = Ref(charge)
    running_sum = Ref(zero(T))

    inactive_layer_exist = !isnothing(ctm.inactive_layer_geometry) || has_inactive_layer(point_types.data)

    bulk_ctmodel = _ctmodel_name(ctm.bulk_charge_trapping_model)
    inactive_layer_ctmodel = _ctmodel_name(ctm.inactive_charge_trapping_model)

    if bulk_ctmodel == "Unknown"
        throw(ArgumentError("Not supported Charge Trapping Model for bulk")) end
    
    if (inactive_layer_ctmodel == "Unknown")||(inactive_layer_ctmodel == "Boggs")
        throw(ArgumentError("Not supported Charge Trapping Model for inactive layer")) end

    τ::T = zero(T)
    τ_inactive::T = zero(T)
    Δt_minimum::T = zero(T)
    
    vth::T = zero(T)
    nσ::T = zero(T)
    
    if bulk_ctmodel == "ConstantLifetime"
        τ = ifelse(charge > 0, ctm.bulk_charge_trapping_model.τh, ctm.bulk_charge_trapping_model.τe)
        Δt_minimum = minimum(diff(pathtimestamps))
        if τ < Δt_minimum
            throw(ArgumentError("The carrier lifetime should be at least bigger than Δt"))
        end
    elseif bulk_ctmodel == "Boggs"
        vth = sqrt(3 * kB * ctm.bulk_charge_trapping_model.temperature / (ifelse(charge > 0, ctm.bulk_charge_trapping_model.meffh, ctm.bulk_charge_trapping_model.meffe) * me)) # in m/s
        nσ = ifelse(charge > 0, ctm.bulk_charge_trapping_model.nσh, ctm.bulk_charge_trapping_model.nσe)
    end

    if inactive_layer_ctmodel == "ConstantLifetime"
        τ_inactive = ifelse(charge > 0, ctm.inactive_charge_trapping_model.τh, ctm.inactive_charge_trapping_model.τe)
        Δt_minimum = minimum(diff(pathtimestamps))
        if τ_inactive < Δt_minimum
            throw(ArgumentError("The carrier lifetime should be at least bigger than Δt"))
        end
    end

    
    @inbounds for i in eachindex(tmp_signal)

        in_inactive_region = in_inactive_layer(path[i], ctm.inactive_layer_geometry, point_types)
        
        Δt::T = (i > 1) ? (pathtimestamps[i] - pathtimestamps[i-1]) : zero(T)
        w::T = i > 1 ? get_interpolation(wpot, path[i], S) : zero(T)

        Δldrift::T = (i > 1) ? norm(path[i] .- path[i-1]) : zero(T)
        Δl::T = Δldrift > 0 ? hypot(Δldrift, vth * (pathtimestamps[i] - pathtimestamps[i-1])) : zero(T)
        
        if in_inactive_region
            if inactive_layer_ctmodel == "ConstantLifetime"
                # Constant trapping lifetime in inactive layer
                tmp_signal[i] = _constantlifetime_signal!( running_sum, τ_inactive, Δt, q, w )
            else
                # No charge trapping in inactive layer
                tmp_signal[i] = _nochargetrapping_signal!(w, charge)
            end

        else
            if bulk_ctmodel == "ConstantLifetime"
                # Constant trapping lifetime in bulk
                tmp_signal[i] = _constantlifetime_signal!( running_sum, τ, Δt, q, w )
            elseif bulk_ctmodel == "Boggs"
                # Boggs trapping model in bulk
                tmp_signal[i] = _Boggs_signal!( running_sum, nσ, Δl, q, w )
            else
                # No charge trapping in bulk
                tmp_signal[i] = _nochargetrapping_signal!(w, charge)
            end
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
    
    bulk_charge_trapping_model = if haskey(config_dict, "model") && config_dict["model"] == "Boggs"
        BoggsChargeTrappingModel{T}(config_dict, temperature = temperature)
    elseif haskey(config_dict, "model") && config_dict["model"] == "ConstantLifetime"
        ConstantLifetimeChargeTrappingModel{T}(config_dict) 
    else
        NoChargeTrappingModel{T}()
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
        NoChargeTrappingModel{T}()
    end
    
    inactive_layer_geometry = haskey(config_dict, "inactive_layer_geometry") ? config_dict["inactive_layer_geometry"] : nothing
    
CombinedChargeTrappingModel{T, typeof(bulk_charge_trapping_model), typeof(inactive_charge_trapping_model), typeof(inactive_layer_geometry)}(
    bulk_charge_trapping_model,
    inactive_charge_trapping_model,
    inactive_layer_geometry
)
end
