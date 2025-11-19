# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).


@doc raw"""
    ADLChargeDriftModel{T <: SSDFloat, M <: AbstractDriftMaterial, N, TM <: AbstractTemperatureModel{T}} <: AbstractChargeDriftModel{T}

Charge drift model for electrons and holes based on the AGATA Detector Library
published in [Bruyneel _et al._, NIMA 569, 764 (2006)](https://doi.org/10.1016/j.nima.2006.08.130).
Find a detailed description of the calculations in [ADL Charge Drift Model](@ref).

## Fields
- `electrons::CarrierParameters{T}`: Parameters to describe the electron drift along the <100> and <111> axes.
- `holes::CarrierParameters{T}`: Parameters to describe the hole drift along the <100> and <111> axes.
- `crystal_orientation::SMatrix{3,3,T,9}`: Rotation matrix that transforms the global coordinate system to the crystal coordinate system given by the <100>, <010> and <001> axes of the crystal.
- `γ::SVector{N,SMatrix{3,3,T,9}}`: Reciprocal mass tensors to the `N` valleys of the conduction band.
- `parameters::ADLParameters{T}`: Parameters needed for the calculation of the electron drift velocity.
- `temperaturemodel::TM`: Models to scale the resulting drift velocities with respect to temperature

See also [`CarrierParameters`](@ref).
"""
struct ADLChargeDriftModel{T <: SSDFloat, M <: AbstractDriftMaterial, N, TM <: AbstractTemperatureModel{T}} <: AbstractChargeDriftModel{T}
    electrons::CarrierParameters{T}
    holes::CarrierParameters{T}
    crystal_orientation::SMatrix{3,3,T,9}
    γ::SVector{N,SMatrix{3,3,T,9}}
    parameters::ADLParameters{T}
    temperaturemodel::TM
end

function _ADLChargeDriftModel(
        config::AbstractDict; 
        T::Type=Float32,
        e100μ0::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["e100"]["mu0"], 
        e100β::Union{RealQuantity,String}  = config["drift"]["velocity"]["parameters"]["e100"]["beta"], 
        e100E0::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["e100"]["E0"],
        e100μn::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["e100"]["mun"],
        e111μ0::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["e111"]["mu0"], 
        e111β::Union{RealQuantity,String}  = config["drift"]["velocity"]["parameters"]["e111"]["beta"], 
        e111E0::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["e111"]["E0"],
        e111μn::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["e111"]["mun"],
        h100μ0::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["h100"]["mu0"], 
        h100β::Union{RealQuantity,String}  = config["drift"]["velocity"]["parameters"]["h100"]["beta"], 
        h100E0::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["h100"]["E0"],
        h111μ0::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["h111"]["mu0"], 
        h111β::Union{RealQuantity,String}  = config["drift"]["velocity"]["parameters"]["h111"]["beta"], 
        h111E0::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["h111"]["E0"],
        input_units::Union{Missing, NamedTuple} = missing, 
        material::Type{<:AbstractDriftMaterial} = HPGe,
        temperature::Union{Missing, Real, Unitful.Temperature} = missing,
        phi110::Union{Missing, Real, AngleQuantity} = missing,
    )::ADLChargeDriftModel{T}
    
    mobility_unit, efield_unit = if ismissing(input_units)
        internal_mobility_unit, internal_efield_unit
    else
        input_units.length^2/(input_units.potential*internal_time_unit), 
        input_units.potential/input_units.length
    end
    # Temperature is not parsed from charge drift config file. It is set by the user or parsed during Semiconductor construction.

    e100 = VelocityParameters{T}(
        _parse_value(T, e100μ0, mobility_unit), 
        _parse_value(T, e100β,  NoUnits), 
        _parse_value(T, e100E0, efield_unit),
        _parse_value(T, e100μn, mobility_unit)
    )
    
    e111 = VelocityParameters{T}(
        _parse_value(T, e111μ0, mobility_unit), 
        _parse_value(T, e111β,  NoUnits), 
        _parse_value(T, e111E0, efield_unit),
        _parse_value(T, e111μn, mobility_unit)
    )
    
    h100 = VelocityParameters{T}(
        _parse_value(T, h100μ0, mobility_unit), 
        _parse_value(T, h100β,  NoUnits), 
        _parse_value(T, h100E0, efield_unit),
        0
    )
    
    h111 = VelocityParameters{T}(
        _parse_value(T, h111μ0, mobility_unit), 
        _parse_value(T, h111β,  NoUnits), 
        _parse_value(T, h111E0, efield_unit),
        0
    )

    electrons = CarrierParameters{T}(e100, e111)
    holes     = CarrierParameters{T}(h100, h111)
    
    
    if "material" in keys(config)
        config_material::String = config["material"]
        if config_material == "HPGe"
            material = HPGe
        elseif config_material == "Si"
            material = Si
        else
            @warn "Material \"$(config_material)\" not supported for ADLChargeDriftModel.\nSupported materials are \"Si\" and \"HPGe\".\nUsing \"$(Symbol(material))\" as default."
        end
    end

    parameters = if "masses" in keys(config["drift"])
        ml = T(config["drift"]["masses"]["ml"])
        mt = T(config["drift"]["masses"]["mt"])
        ADLParameters{T}(ml, mt, material)
    else
        ml = T(material_properties[Symbol(material)].ml)
        mt = T(material_properties[Symbol(material)].mt)
        ADLParameters{T}(ml, mt, material)
    end

    crystal_orientation::SMatrix{3,3,T,9} = if ismissing(phi110) 
        if "phi110" in keys(config)
            RotZ{T}(- π/4 - _parse_value(T, config["phi110"], u"rad"))
        else
            transpose(parse_rotation_matrix(T, config, u"rad")) # replace u"rad" with the units in the config file
        end
    else
        RotZ{T}(- π/4 - _parse_value(T, phi110, u"rad"))
    end

    γ = setup_γj(crystal_orientation, parameters, material)                      

    temperaturemodel = TemperatureModel(T, config)

    cdm = ADLChargeDriftModel{T,material,length(γ),typeof(temperaturemodel)}(electrons, holes, crystal_orientation, γ, parameters, temperaturemodel)

    scale_to_temperature(cdm, temperature)
end

# Check the syntax of the ADLChargeDriftModel config file before parsing
function ADLChargeDriftModel(config::AbstractDict, input_units::Union{Missing, NamedTuple} = missing; kwargs...)
    if !haskey(config, "drift") 
        throw(ConfigFileError("ADLChargeDriftModel config file needs entry 'drift'.")) 
    elseif !haskey(config["drift"], "velocity") 
        throw(ConfigFileError("ADLChargeDriftModel config file needs entry 'drift/velocity'.")) 
    elseif !haskey(config["drift"]["velocity"], "parameters")
        throw(ConfigFileError("ADLChargeDriftModel config file needs entry 'drift/velocity/parameters'."))
    end

    for axis in ("e100", "e111", "h100", "h111")
        if !haskey(config["drift"]["velocity"]["parameters"], axis)
            throw(ConfigFileError("ADLChargeDriftModel config file needs entry 'drift/velocity/parameters/$(axis)'."))
        end
        for param in ("mu0", "beta", "E0", "mun")
            if !haskey(config["drift"]["velocity"]["parameters"][axis], param) && !(axis[1] == 'h' && param == "mun") # holes have no μn
                throw(ConfigFileError("ADLChargeDriftModel config file needs entry 'drift/velocity/parameters/$(axis)/$(param)'."))
            end
        end
    end
    _ADLChargeDriftModel(config; input_units = input_units, kwargs...)
end

const default_ADL_config_file = joinpath(get_path_to_example_config_files(), "ADLChargeDriftModel/drift_velocity_config.yaml")
function ADLChargeDriftModel(config_filename::AbstractString = default_ADL_config_file, input_units::Union{Missing, NamedTuple} = missing; kwargs...)
    ADLChargeDriftModel(parse_config_file(config_filename), input_units; kwargs...)
end

ADLChargeDriftModel{T}(args...; kwargs...) where {T <: SSDFloat} = ADLChargeDriftModel(args...; T=T, kwargs...)