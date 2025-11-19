"""
    struct IsotropicChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
        
Charge drift model in which the electrons and holes drift along the electric field with a constant mobility in m²/Vs

## Fields
- `μ_e::Union{RealQuantity, AbstractString}`: Mobility of the electrons in m²/Vs.
- `μ_h::Union{RealQuantity, AbstractString}`: Mobility of the holes in m²/Vs.

"""
struct IsotropicChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
    μ_e::T
    μ_h::T
end

function _IsotropicChargeDriftModel(
        config::AbstractDict; 
        T::Type=Float32,
        μ_e::Union{RealQuantity, AbstractString} = config["mobilities"]["e"], 
        μ_h::Union{RealQuantity, AbstractString} = config["mobilities"]["h"],
        input_units::Union{Missing, NamedTuple} = missing,
        temperature::Union{Missing, Real, Unitful.Temperature} = missing
    )::IsotropicChargeDriftModel{T}
    mobility_unit = ismissing(input_units) ? internal_mobility_unit : input_units.length^2/(input_units.potential*internal_time_unit)
   IsotropicChargeDriftModel{T}(_parse_value.(T, (μ_e, μ_h), mobility_unit)...)
end


# Check the syntax of the IsotropicChargeDriftModel config file before parsing
function IsotropicChargeDriftModel(config::AbstractDict, input_units::Union{Missing, NamedTuple} = missing; kwargs...)
    if !haskey(config, "mobilities")
        throw(ConfigFileError("IsotropicChargeDriftModel config file needs entry 'mobilities'."))
    end
    for axis in ("e", "h")
        if !haskey(config["mobilities"], axis)
            throw(ConfigFileError("IsotropicChargeDriftModel config file needs entry 'mobilities/$axis'."))
        end
    end
    _IsotropicChargeDriftModel(config; input_units = input_units, kwargs...)
end

const default_icd_config_file = joinpath(get_path_to_example_config_files(), "IsotropicChargeDriftModel/mobilities.yaml")
function IsotropicChargeDriftModel(config_filename::AbstractString = default_icd_config_file, input_units::Union{Missing, NamedTuple} = missing; kwargs...)
    IsotropicChargeDriftModel(parse_config_file(config_filename), input_units; kwargs...)
end

IsotropicChargeDriftModel{T}(args...; kwargs...) where {T <: SSDFloat} = IsotropicChargeDriftModel(args...; T=T, kwargs...)

calculate_mobility(cdm::IsotropicChargeDriftModel{T}, ::AbstractCoordinatePoint{T}, ::Type{Electron}) where {T <: SSDFloat} = cdm.μ_e
calculate_mobility(cdm::IsotropicChargeDriftModel{T}, ::AbstractCoordinatePoint{T}, ::Type{Hole})     where {T <: SSDFloat} = cdm.μ_h

@fastmath function getVe(fv::SVector{3, T}, cdm::IsotropicChargeDriftModel{T}, current_pos::CartesianPoint{T} = zero(CartesianPoint{T})) where {T <: SSDFloat}
    -cdm.μ_e*fv
end

@fastmath function getVh(fv::SVector{3, T}, cdm::IsotropicChargeDriftModel{T}, current_pos::CartesianPoint{T} = zero(CartesianPoint{T})) where {T <: SSDFloat}
    cdm.μ_h*fv
end