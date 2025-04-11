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
    IsotropicChargeDriftModel{T}(
        μ_e::Union{RealQuantity, AbstractString}, 
        μ_h::Union{RealQuantity, AbstractString}
    ) where {T <: SSDFloat} = new{T}(_parse_value.(T, (μ_e, μ_h), internal_mobility_unit)...)
end

IsotropicChargeDriftModel(args...; T::Type = Float32, kwargs...) = IsotropicChargeDriftModel{T}(args...; kwargs...)

const default_icd_config_file = joinpath(get_path_to_example_config_files(), "IsotropicChargeDriftModel/mobilities.yaml")

function IsotropicChargeDriftModel{T}(config_filename::AbstractString = default_icd_config_file) where {T <: SSDFloat}
    IsotropicChargeDriftModel(parse_config_file(config_filename))
end

function IsotropicChargeDriftModel{T}(config::AbstractDict) where {T <: SSDFloat}
    if !haskey(config, "mobilities")
        throw(ConfigFileError("IsotropicChargeDriftModel config file needs entry 'mobilities'."))
    end

    for axis in ("e", "h")
        if !haskey(config["mobilities"], axis)
            throw(ConfigFileError("IsotropicChargeDriftModel config file needs entry 'mobilities/$axis'."))
        end
    end
    
    IsotropicChargeDriftModel{T}(config["mobilities"]["e"], config["mobilities"]["h"])
end

@fastmath function getVe(fv::SVector{3, T}, cdm::IsotropicChargeDriftModel{T}) where {T <: SSDFloat}
    @inbounds begin 
        -cdm.μ_e*fv
    end
end

@fastmath function getVh(fv::SVector{3, T}, cdm::IsotropicChargeDriftModel{T}) where {T <: SSDFloat}
    @inbounds begin 
        cdm.μ_h*fv
    end
end