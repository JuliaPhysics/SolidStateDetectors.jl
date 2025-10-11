# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

@doc raw"""
    ADL2016ChargeDriftModel{T <: SSDFloat, M <: AbstractDriftMaterial, N, TM <: AbstractTemperatureModel{T}} <: AbstractChargeDriftModel{T}

Charge drift model for electrons and holes based on the AGATA Detector Library,
published in [Bruyneel _et al._, EPJA 52, 70 (2016)](https://doi.org/10.1140/epja/i2016-16070-9).

## Fields
- `electrons::VelocityParameters{T}`: Parameters to describe the electron drift along the <100> axis.
- `holes::CarrierParameters{T}`: Parameters to describe the hole drift along the <100> and <111> axes.
- `crystal_orientation::SMatrix{3,3,T,9}`: Rotation matrix that transforms the global coordinate system to the crystal coordinate system given by the <100>, <010> and <001> axes of the crystal.
- `γ::SVector{N,SMatrix{3,3,T,9}}`: Reciprocal mass tensors to the `N` valleys of the conduction band.
- `parameters::ADLParameters{T}`: Parameters needed for the calculation of the electron drift velocity.
- `temperaturemodel::TM`: Models to scale the resulting drift velocities with respect to temperature

See also [`VelocityParameters`](@ref) and [`CarrierParameters`](@ref).
"""
struct ADL2016ChargeDriftModel{T <: SSDFloat, M <: AbstractDriftMaterial, N, TM <: AbstractTemperatureModel{T}} <: AbstractChargeDriftModel{T}
    electrons::VelocityParameters{T}
    holes::CarrierParameters{T}
    crystal_orientation::SMatrix{3,3,T,9}
    γ::SVector{N,SMatrix{3,3,T,9}}
    parameters::ADLParameters{T}
    temperaturemodel::TM
end

function _ADL2016ChargeDriftModel(
        config::AbstractDict; 
        T::Type=Float32,
        e100μ0::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["e100"]["mu0"], 
        e100β::Union{RealQuantity,String}  = config["drift"]["velocity"]["parameters"]["e100"]["beta"], 
        e100E0::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["e100"]["E0"],
        e100μn::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["e100"]["mun"],
        η0::Union{Real,String}             = config["drift"]["velocity"]["parameters"]["escattering"]["eta0"],
        b::Union{Real,String}              = config["drift"]["velocity"]["parameters"]["escattering"]["b"],
        Eref::Union{RealQuantity,String}   = config["drift"]["velocity"]["parameters"]["escattering"]["Eref"],
        h100μ0::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["h100"]["mu0"], 
        h100β::Union{RealQuantity,String}  = config["drift"]["velocity"]["parameters"]["h100"]["beta"], 
        h100E0::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["h100"]["E0"],
        h111μ0::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["h111"]["mu0"], 
        h111β::Union{RealQuantity,String}  = config["drift"]["velocity"]["parameters"]["h111"]["beta"], 
        h111E0::Union{RealQuantity,String} = config["drift"]["velocity"]["parameters"]["h111"]["E0"],
        material::Type{<:AbstractDriftMaterial} = HPGe,
        temperature::Union{Missing, Real} = missing, 
        phi110::Union{Missing, Real, AngleQuantity} = missing
    )::ADL2016ChargeDriftModel{T}
    
    electrons = VelocityParameters{T}(
        _parse_value(T, e100μ0, internal_mobility_unit), 
        _parse_value(T, e100β,  NoUnits), 
        _parse_value(T, e100E0, internal_efield_unit),
        _parse_value(T, e100μn, internal_mobility_unit)
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

    ml_inv, mt_inv = parameters = if "masses" in keys(config["drift"])
        ml = T(config["drift"]["masses"]["ml"])
        mt = T(config["drift"]["masses"]["mt"])
        inv(ml), inv(mt)
    else
        ml = T(material_properties[Symbol(material)].ml)
        mt = T(material_properties[Symbol(material)].mt)
        inv(ml), inv(mt)
    end
    
    parameters = ADLParameters{T}(
        ml_inv, mt_inv, 
        _parse_value(T, η0, NoUnits), 
        _parse_value(T, b, NoUnits), 
        _parse_value(T, Eref, internal_efield_unit)
    )

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
    return ADL2016ChargeDriftModel{T,material,length(γ),typeof(temperaturemodel)}(electrons, holes, crystal_orientation, γ, parameters, temperaturemodel)
end

# Check the syntax of the ADL2016ChargeDriftModel config file before parsing
function ADL2016ChargeDriftModel(config::AbstractDict; kwargs...)
    if !haskey(config, "drift") 
        throw(ConfigFileError("ADL2016ChargeDriftModel config file needs entry 'drift'.")) 
    elseif !haskey(config["drift"], "velocity") 
        throw(ConfigFileError("ADL2016ChargeDriftModel config file needs entry 'drift/velocity'.")) 
    elseif !haskey(config["drift"]["velocity"], "parameters")
        throw(ConfigFileError("ADL2016ChargeDriftModel config file needs entry 'drift/velocity/parameters'."))
    end

    for axis in ("e100", "escattering", "h100", "h111")
        if !haskey(config["drift"]["velocity"]["parameters"], axis)
            throw(ConfigFileError("ADLCharge2016DriftModel config file needs entry 'drift/velocity/parameters/$(axis)'."))
        end
        for param in (axis == "escattering" ? ("eta0", "b", "Eref") : ("mu0", "beta", "E0", "mun"))
            if !haskey(config["drift"]["velocity"]["parameters"][axis], param) && !(axis[1] == 'h' && param == "mun") # holes have no μn
                throw(ConfigFileError("ADLCharge2016DriftModel config file needs entry 'drift/velocity/parameters/$(axis)/$(param)'."))
            end
        end
    end
    _ADL2016ChargeDriftModel(config; kwargs...)
end

const default_ADL2016_config_file = joinpath(get_path_to_example_config_files(), "ADLChargeDriftModel/drift_velocity_config_2016.yaml")
function ADL2016ChargeDriftModel(config_filename::AbstractString = default_ADL2016_config_file; kwargs...)
    ADL2016ChargeDriftModel(parse_config_file(config_filename); kwargs...)
end

ADL2016ChargeDriftModel{T}(args...; kwargs...) where {T <: SSDFloat} = ADL2016ChargeDriftModel(args...; T=T, kwargs...)
