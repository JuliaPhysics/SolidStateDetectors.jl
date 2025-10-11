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

function ADLChargeDriftModel{T,M,N,TM}(chargedriftmodel::ADLChargeDriftModel{<:Any,M,N,TM})::ADLChargeDriftModel{T,M,N,TM} where {T <: SSDFloat, M, N, TM}
    cdmf64 = chargedriftmodel
    e100 = VelocityParameters{T}(cdmf64.electrons.axis100.mu0, cdmf64.electrons.axis100.beta, cdmf64.electrons.axis100.E0, cdmf64.electrons.axis100.mun)
    e111 = VelocityParameters{T}(cdmf64.electrons.axis111.mu0, cdmf64.electrons.axis111.beta, cdmf64.electrons.axis111.E0, cdmf64.electrons.axis111.mun)
    h100 = VelocityParameters{T}(cdmf64.holes.axis100.mu0, cdmf64.holes.axis100.beta, cdmf64.holes.axis100.E0, cdmf64.holes.axis100.mun)
    h111 = VelocityParameters{T}(cdmf64.holes.axis111.mu0, cdmf64.holes.axis111.beta, cdmf64.holes.axis111.E0, cdmf64.holes.axis111.mun)
    electrons = CarrierParameters{T}(e100, e111)
    holes     = CarrierParameters{T}(h100, h111)
    crystal_orientation = SMatrix{3,3,T,9}(cdmf64.crystal_orientation)
    γ = Vector{SMatrix{3,3,T,9}}( cdmf64.γ )
    parameters = ADLParameters{T}(cdmf64.parameters.ml_inv, cdmf64.parameters.mt_inv, cdmf64.parameters.Γ0_inv, cdmf64.parameters.Γ1, cdmf64.parameters.Γ2)
    temperaturemodel::AbstractTemperatureModel{T} = cdmf64.temperaturemodel
    ADLChargeDriftModel{T,M,N,TM}(electrons, holes, crystal_orientation, γ, parameters, temperaturemodel)
end



ADLChargeDriftModel{T}(args...; kwargs...) where {T <: SSDFloat} = ADLChargeDriftModel(args...; T=T, kwargs...)

const default_ADL_config_file = joinpath(get_path_to_example_config_files(), "ADLChargeDriftModel/drift_velocity_config.yaml")
function ADLChargeDriftModel(config_filename::AbstractString = default_ADL_config_file; kwargs...)
    ADLChargeDriftModel(parse_config_file(config_filename); kwargs...)
end

# Check the syntax of the ADLChargeDriftModel config file before parsing
function ADLChargeDriftModel(config::AbstractDict; kwargs...)
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
    _ADLChargeDriftModel(config; kwargs...)
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
        material::Type{<:AbstractDriftMaterial} = HPGe,
        temperature::Union{Missing, Real}= missing, 
        phi110::Union{Missing, Real, AngleQuantity} = missing
    )::ADLChargeDriftModel{T}
    
    e100 = VelocityParameters{T}(
        _parse_value(T, e100μ0, internal_mobility_unit), 
        _parse_value(T, e100β,  NoUnits), 
        _parse_value(T, e100E0, internal_efield_unit),
        _parse_value(T, e100μn, internal_mobility_unit)
    )
    
    e111 = VelocityParameters{T}(
        _parse_value(T, e111μ0, internal_mobility_unit), 
        _parse_value(T, e111β,  NoUnits), 
        _parse_value(T, e111E0, internal_efield_unit),
        _parse_value(T, e111μn, internal_mobility_unit)
    )
    
    h100 = VelocityParameters{T}(
        _parse_value(T, h100μ0, internal_mobility_unit), 
        _parse_value(T, h100β,  NoUnits), 
        _parse_value(T, h100E0, internal_efield_unit),
        0
    )
    
    h111 = VelocityParameters{T}(
        _parse_value(T, h111μ0, internal_mobility_unit), 
        _parse_value(T, h111β,  NoUnits), 
        _parse_value(T, h111E0, internal_efield_unit),
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

    if !ismissing(temperature) temperature = T(temperature) end  #if you give the temperature it will be used, otherwise read from config file                         

    if "temperature_dependence" in keys(config)
        if "model" in keys(config["temperature_dependence"])
            model::String = config["temperature_dependence"]["model"]
            if model == "Linear"
                temperaturemodel = LinearModel{T}(config, temperature = temperature)
            elseif model == "PowerLaw"
                temperaturemodel = PowerLawModel{T}(config, temperature = temperature)
            elseif model == "Boltzmann"
                temperaturemodel = BoltzmannModel{T}(config, temperature = temperature)
            else
                temperaturemodel = VacuumModel{T}(config)
                println("Config File does not suit any of the predefined temperature models. The drift velocity will not be rescaled.")
            end
        else
            temperaturemodel = VacuumModel{T}(config)
            println("No temperature model specified. The drift velocity will not be rescaled.")
        end
    else
        temperaturemodel = VacuumModel{T}(config)
        # println("No temperature dependence found in Config File. The drift velocity will not be rescaled.")
    end
    return ADLChargeDriftModel{T,material,length(γ),typeof(temperaturemodel)}(electrons, holes, crystal_orientation, γ, parameters, temperaturemodel)
end
