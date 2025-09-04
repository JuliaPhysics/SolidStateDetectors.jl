abstract type AbstractChargeTrappingModel{T <: SSDFloat} end

function _calculate_signal( 
        ::AbstractChargeTrappingModel{T},
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
    struct ConstantLifetimeChargeTrappingModelInactiveLayer{T <: SSDFloat} <: AbstractChargeTrappingModel{T}
        
This constant-lifetime-based charge trapping model is similar to the Boggs model, which is constant-mean-free-path based.

## Fields
* `τh_inactive::T`: Lifetime for holes in surface layer (dead layer/inactive layer) (default: `τh_inactive = 80ns`).
* `τe_inactive::T`: Lifetime for electrons in surface layer (default: `τe_inactive = 80ns`).
* `inactive_layer_geometry`: The geometry of the inactive layer. If it's not defined, it uses point_types to get the inactive layer geometry, else it uses the user-defined inactive_layer_geometry.
> Note: All τ must be much bigger than `Δt`.

See also [Charge Trapping Models](@ref).
"""

struct ConstantLifetimeChargeTrappingModelInactiveLayer{T <: SSDFloat, G <: Union{<:AbstractGeometry, Nothing}} <: AbstractChargeTrappingModel{T}
    τh_inactive::T
    τe_inactive::T
    inactive_layer_geometry::G
end

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

ConstantLifetimeChargeTrappingModelInactiveLayer(args...; T::Type{<:SSDFloat}, kwargs...) = ConstantLifetimeChargeTrappingModelInactiveLayer{T}(args...; kwargs...)
function ConstantLifetimeChargeTrappingModelInactiveLayer{T}(config_dict::AbstractDict = Dict()) where {T <: SSDFloat}
    τh_inactive::T = ustrip(u"s", 80u"ns")
    τe_inactive::T = ustrip(u"s", 80u"ns")

    if haskey(config_dict, "model_inactive") && !haskey(config_dict, "parameters_inactive")
        throw(ConfigFileError("`ConstantLifetimeChargeTrappingModelInactiveLayer` does not have `parameters_inactive`"))
    end

    parameters_inactive = haskey(config_dict, "parameters_inactive") ? config_dict["parameters_inactive"] : config_dict
    inactive_layer_geometry = haskey(config_dict, "inactive_layer_geometry") ? config_dict["inactive_layer_geometry"] : nothing

    allowed_keys = ("τh_inactive", "τe_inactive")
    k = filter(k -> !(k in allowed_keys), keys(parameters_inactive))
    !isempty(k) && @warn "The following keys will be ignored: $(k).\nAllowed keys are: $(allowed_keys)"

    if haskey(parameters_inactive, "τh_inactive")   τh_inactive =     _parse_value(T, parameters_inactive["τh_inactive"], internal_time_unit) end
    if haskey(parameters_inactive, "τe_inactive")   τe_inactive =     _parse_value(T, parameters_inactive["τe_inactive"], internal_time_unit) end
    ConstantLifetimeChargeTrappingModelInactiveLayer{T, typeof(inactive_layer_geometry)}(τh_inactive, τe_inactive, inactive_layer_geometry)
end


"""
    struct BoggsChargeTrappingModel{T <: SSDFloat} <: AbstractChargeTrappingModel{T}
        
Charge trapping model presented in [Boggs _et al._ (2023)](https://doi.org/10.1016/j.nima.2023.168756).

## Fields
* `nσe::T`: Trapping product for electrons (default: `(nσe)^-1 = 1020cm`).
* `nσh::T`: Trapping product for holes (default: `(nσh)^-1 = 2040cm`).
* `temperature::T`: Temperature of the crystal (default: `78K`).
* `inactive_layer_geometry`: The geometry of the inactive layer. If it's not defined, it uses point_types to get the inactive layer geometry, else it uses the user-defined inactive_layer_geometry.

See also [Charge Trapping Models](@ref).
"""
struct BoggsChargeTrappingModel{T <: SSDFloat, G <: Union{<:AbstractGeometry, Nothing}} <: AbstractChargeTrappingModel{T}
    nσe::T  # in m^-1
    nσh::T  # in m^-1
    meffe::T # in units of me
    meffh::T # in units of me
    temperature::T # in K
    inactive_layer_geometry::G
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
        ctmil::AbstractChargeTrappingModel{T},
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

    inactive_cct_defined = (ctmil isa ConstantLifetimeChargeTrappingModelInactiveLayer) &&
            hasproperty(ctmil, :τh_inactive) &&
            hasproperty(ctmil, :τe_inactive)

    vth::T = sqrt(3 * kB * ctm.temperature / (ifelse(charge > 0, ctm.meffh, ctm.meffe) * me)) # in m/s
    nσ::T = ifelse(charge > 0, ctm.nσh, ctm.nσe)
    Δt_minimum = minimum(diff(pathtimestamps))
    τ_inactive = if inactive_cct_defined
        ifelse(charge > 0, ctmil.τh_inactive, ctmil.τe_inactive)::T
    else
        nothing
    end

    if (inactive_cct_defined && τ_inactive<Δt_minimum)
        throw(ArgumentError("The carrier lifetime should be at least bigger than Δt"))
    end

    @inbounds for i in eachindex(tmp_signal)

        in_inactive_region = in_inactive_layer(path[i], ctm.inactive_layer_geometry, point_types)
        
        Δt::T = (i > 1) ? (pathtimestamps[i] - pathtimestamps[i-1]) : zero(T)
        Δldrift::T = (i > 1) ? norm(path[i] .- path[i-1]) : zero(T)
        Δl::T = Δldrift > 0 ? hypot(Δldrift, vth * (pathtimestamps[i] - pathtimestamps[i-1])) : zero(T)
        w::T = i > 1 ? get_interpolation(wpot, path[i], S) : zero(T)

        
        if in_inactive_region && inactive_cct_defined
            # Constant trapping lifetime in inactive layer
            tmp_signal[i] = _constantlifetime_signal!( running_sum, τ_inactive, Δt, q, w )

        elseif in_inactive_region && (typeof(ctmil) <: NoChargeTrappingModel)
            # No charge trapping in inactive layer
            tmp_signal[i] = _nochargetrapping_signal!(w, charge)

        elseif in_inactive_region
            throw("For the chosen charge trapping model in the bulk, no `_calculate_signal` is implemented for the chosen charge trapping model in the inactive layer.")

        else
            # Boggs trapping model in bulk
            tmp_signal[i] = _Boggs_signal!( running_sum, nσ, Δl, q, w )

        end
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

    inactive_layer_geometry = haskey(config_dict, "inactive_layer_geometry") ? config_dict["inactive_layer_geometry"] : nothing
    
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
    BoggsChargeTrappingModel{T, typeof(inactive_layer_geometry)}(nσe, nσh, meffe, meffh, temperature, inactive_layer_geometry)
end




"""
    struct ConstantLifetimeChargeTrappingModel{T <: SSDFloat} <: AbstractChargeTrappingModel{T}
        
This constant-lifetime-based charge trapping model is similar to the Boggs model, which is constant-mean-free-path based.

## Fields
* `τh::T`: Lifetime for holes in bulk volume (sensitive region) (default: `τh = 1ms`).
* `τe::T`: Lifetime for electrons in bulk volume (default: `τe = 1ms`).
* `inactive_layer_geometry`: The geometry of the inactive layer. If it's not defined, it uses point_types to get the inactive layer geometry, else it uses the user-defined inactive_layer_geometry.
> Note: All τ must be much bigger than `Δt`.

See also [Charge Trapping Models](@ref).
"""

struct ConstantLifetimeChargeTrappingModel{T <: SSDFloat, G <: Union{<:AbstractGeometry, Nothing}} <: AbstractChargeTrappingModel{T}
    τh::T
    τe::T
    inactive_layer_geometry::G
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
        ctmil::AbstractChargeTrappingModel{T},
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

    inactive_cct_defined = (ctmil isa ConstantLifetimeChargeTrappingModelInactiveLayer) &&
            hasproperty(ctmil, :τh_inactive) &&
            hasproperty(ctmil, :τe_inactive)


    τ::T = ifelse(charge > 0, ctm.τh, ctm.τe)
    Δt_minimum = minimum(diff(pathtimestamps))
    τ_inactive = if inactive_cct_defined
        ifelse(charge > 0, ctmil.τh_inactive, ctmil.τe_inactive)::T
    else
        nothing
    end

    if τ<Δt_minimum || (inactive_cct_defined && τ_inactive<Δt_minimum) 
        throw(ArgumentError("The carrier lifetime should be at least bigger than Δt"))
    end
   
    
    @inbounds for i in eachindex(tmp_signal)

        in_inactive_region = in_inactive_layer(path[i], ctm.inactive_layer_geometry, point_types)
        
        Δt::T = (i > 1) ? (pathtimestamps[i] - pathtimestamps[i-1]) : zero(T)
        w::T = i > 1 ? get_interpolation(wpot, path[i], S) : zero(T)

        if in_inactive_region && inactive_cct_defined
            # Constant trapping lifetime in inactive layer
            tmp_signal[i] = _constantlifetime_signal!( running_sum, τ_inactive, Δt, q, w )

        elseif in_inactive_region && (typeof(ctmil) <: NoChargeTrappingModel)
            # No charge trapping in inactive layer
            tmp_signal[i] = _nochargetrapping_signal!(w, charge)

        elseif in_inactive_region
            throw("For the chosen charge trapping model in the bulk, no `_calculate_signal` is implemented for the chosen charge trapping model in the inactive layer.")

        else
            # Constant trapping lifetime in bulk
            tmp_signal[i] = _constantlifetime_signal!( running_sum, τ, Δt, q, w )

        end
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
    inactive_layer_geometry = haskey(config_dict, "inactive_layer_geometry") ? config_dict["inactive_layer_geometry"] : nothing

    allowed_keys = ("τh", "τe")
    k = filter(k -> !(k in allowed_keys), keys(parameters))
    !isempty(k) && @warn "The following keys will be ignored: $(k).\nAllowed keys are: $(allowed_keys)"

    if haskey(parameters, "τh")   τh =     _parse_value(T, parameters["τh"], internal_time_unit) end
    if haskey(parameters, "τe")   τe =     _parse_value(T, parameters["τe"], internal_time_unit) end
    ConstantLifetimeChargeTrappingModel{T, typeof(inactive_layer_geometry)}(τh, τe, inactive_layer_geometry)
end
