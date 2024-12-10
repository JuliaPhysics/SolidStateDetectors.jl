@doc raw"""

    struct VelocityParameters{T <: SSDFloat}

Values needed to parametrize the longitudinal drift velocity of electrons or hole
along a crystal axis as a function of the electric field strength.

## Background information 

The parameterization for the longitudinal drift velocity, ``v_l``, as a function of the electric 
field strength, ``E``, was proposed by [D.M. Caughey and R.E. Thomas](https://ieeexplore.ieee.org/document/1448053)
and later expanded by [L. Mihailescu et al.](https://www.sciencedirect.com/science/article/pii/S0168900299012863):
```math
v_l = \frac{\mu_0 E}{(1 + (E/E_0 )^{\beta})^{1/ \beta}} - \mu_{n} E.
```
with the four parameters, ``\mu_0``, ``\beta``, ``E_0`` and ``\mu_n``, which are different
for electrons and holes and for the different crystal axes.

!!! note
    The parameter ``\mu_n`` accounts for the Gunn effects for electrons and should be 0 for holes.
    
## Fields
* `mu0::T`: Parameter ``\mu_0`` in the parameterization shown above.
* `beta::T`: Parameter ``\beta`` in the parameterization shown above.
* `E0::T`: Parameter ``E_0`` in the parameterization shown above.
* `mun::T`: Parameter ``\mu_n`` in the parameterization shown above.

"""
struct VelocityParameters{T <: SSDFloat}
    mu0::T
    beta::T
    E0::T
    mun::T
end

@doc raw"""

    struct CarrierParameters{T <: SSDFloat}

Parameters needed to describe the electron or hole drift along the <100> and <111> axes.

## Fields
* `axis100::VelocityParameters{T}`: Parameters to describe the charge drift along the <100> axis.
* `axis111::VelocityParameters{T}`: Parameters to describe the charge drift along the <111> axis.

See also [`VelocityParameters`](@ref).
"""
struct CarrierParameters{T <: SSDFloat}
    axis100::VelocityParameters{T}
    axis111::VelocityParameters{T}
end

struct ADLParameters{T <: SSDFloat}
    ml_inv::T
    mt_inv::T
    Γ0::T
    Γ1::T
    Γ2::T
end

# Temperature models
include("TemperatureModels/TemperatureModels.jl")

# Electron model parametrization from [3]
@fastmath function γj(j::Integer, crystal_orientation::SMatrix{3,3,T,9}, γ0::SMatrix{3,3,T,9}, ::Type{HPGe})::SMatrix{3,3,T,9} where {T <: SSDFloat}
    tmp::T = 2 / 3
    a::T = acos(sqrt(tmp))
    Rx::SMatrix{3,3,T,9} = RotX{T}(a)
    b::T = (j - 1) * T(π) / 2 - π / 4
    Rzj::SMatrix{3,3,T,9} = RotZ{T}(-b)
    Rj = Rx * Rzj * crystal_orientation
    transpose(Rj) * γ0 * Rj
end

@fastmath function γj(j::Integer, crystal_orientation::SMatrix{3,3,T,9}, γ0::SMatrix{3,3,T,9}, ::Type{Si})::SMatrix{3,3,T,9} where {T <: SSDFloat}
    # needs to be updated!
    b::T = (j == 1 ? 0 : π / 2)
    Rzj::SMatrix{3,3,T,9} = (j == 3 ? SMatrix{3,3,T,9}(1, 0, 0, 0, 0, 1, 0, -1, 0)  : I) * SMatrix{3,3,T,9}(cos(b), -sin(b), 0, sin(b), cos(b), 0, 0, 0, 1) 
    Rj = Rzj * crystal_orientation
    transpose(Rj) * γ0 * Rj
end

@fastmath function setup_γj(crystal_orientation::SMatrix{3,3,T,9}, p::ADLParameters{T}, material::Type{HPGe})::SVector{4, SMatrix{3,3,T,9}} where {T <: SSDFloat}
    γ0 = SMatrix{3,3,T,9}(p.mt_inv, 0, 0, 0, p.ml_inv, 0, 0, 0, p.mt_inv)
    SVector{4, SMatrix{3,3,T,9}}([γj(i, crystal_orientation, γ0, material) for i in Base.OneTo(4)]...)
end

@fastmath function setup_γj(crystal_orientation::SMatrix{3,3,T,9}, p::ADLParameters{T}, material::Type{Si})::SVector{3, SMatrix{3,3,T,9}} where {T <: SSDFloat}
    γ0 = SMatrix{3,3,T,9}(p.mt_inv, 0, 0, 0, p.ml_inv, 0, 0, 0, p.mt_inv)
    SVector{3, SMatrix{3,3,T,9}}([γj(i, crystal_orientation, γ0, material) for i in Base.OneTo(3)]...)
end

@fastmath function ADLParameters{T}(ml::T, mt::T, ::Type{HPGe})::ADLParameters{T} where {T}
    ml_inv::T = inv(ml)
    mt_inv::T = inv(mt)
    tmp::T = sqrt(8*mt_inv + ml_inv) / 3
    
    Γ0::T = inv(sqrt((2 * mt_inv + ml_inv) / 3))
    Γ1::T = (-4 * sqrt(ml_inv) - 4/3 * tmp) / (sqrt(ml_inv) - tmp)^2
    Γ2::T = Γ1 * (sqrt(ml_inv) + 3 * tmp) / (-4)
    ADLParameters{T}(ml_inv, mt_inv, Γ0, Γ1, Γ2)
end

@fastmath function ADLParameters{T}(ml::T, mt::T, ::Type{Si})::ADLParameters{T} where {T}
    ml_inv::T = inv(ml)
    mt_inv::T = inv(mt)
    
    Γ0::T = inv(sqrt((2 * mt_inv + ml_inv) / 3))
    Γ1::T = (-3 * sqrt(ml_inv) - 3/2 * sqrt(mt_inv) )/(sqrt(ml_inv) - sqrt(mt_inv))^2
    Γ2::T = Γ1 * (sqrt(ml_inv) + 2 * sqrt(mt_inv)) / (-3)
    ADLParameters{T}(ml_inv, mt_inv, Γ0, Γ1, Γ2)
end


# Longitudinal drift velocity formula
@fastmath function Vl(Emag::T, params::VelocityParameters{T})::T where {T <: SSDFloat}
    params.mu0 * Emag / (1 + (Emag / params.E0)^params.beta)^(1 / params.beta) - params.mun * Emag
end


@doc raw"""
    ADLChargeDriftModel{T <: SSDFloat, M <: AbstractDriftMaterial, N, TM <: AbstractTemperatureModel{T}} <: AbstractChargeDriftModel{T}

Charge drift model for electrons and holes based on the AGATA Detector Library.
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
            elseif model == "SquareRoot"
                temperaturemodel = SquareRootModel{T}(config, temperature = temperature)
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


@inline function _get_AE_and_RE(cdm::ADLChargeDriftModel{T, HPGe}, V100e::T, V111e::T)::Tuple{T,T} where{T <: SSDFloat}
    edmp::ADLParameters{T} = cdm.parameters
    AE::T = edmp.Γ0 * V100e 
    RE::T = edmp.Γ1 * V111e / AE + edmp.Γ2
    AE, RE
end

@inline function _get_AE_and_RE(cdm::ADLChargeDriftModel{T, Si}, V100e::T, V111e::T)::Tuple{T,T} where{T <: SSDFloat}
    edmp::ADLParameters{T} = cdm.parameters
    AE::T = edmp.Γ0 * V111e
    RE::T = edmp.Γ1 * V100e / AE + edmp.Γ2
    AE, RE
end


@fastmath function getVe(fv::SVector{3, T}, cdm::ADLChargeDriftModel{T,M,N}, Emag_threshold::T = T(1e-5))::SVector{3, T} where {T <: SSDFloat, M, N}
    @inbounds begin
        Emag::T = norm(fv)
        Emag_inv::T = inv(Emag)

        if Emag < Emag_threshold return SVector{3,T}(0, 0, 0) end

        f::NTuple{4,T} = scale_to_given_temperature(Emag, cdm.temperaturemodel)
        f100e::T = f[1]
        f111e::T = f[2]
        V100e::T = Vl(Emag, cdm.electrons.axis100) * f100e
        V111e::T = Vl(Emag, cdm.electrons.axis111) * f111e
        
        AE::T, RE::T = _get_AE_and_RE(cdm, V100e, V111e)

        E0 = SVector{3, T}(fv * Emag_inv)
        oneOverSqrtEγE::SVector{N,T} = broadcast(γ -> T(1/sqrt(γ * E0 ⋅ E0)), cdm.γ)
        sumOneOverSqrtEγE_inv::T = inv(sum(oneOverSqrtEγE))
        N_inv::T = T(1/N)

        g0::SVector{3,T} = @SVector T[0,0,0]
        for j in eachindex(cdm.γ)
            NiOverNj::T = RE * (oneOverSqrtEγE[j] * sumOneOverSqrtEγE_inv - N_inv) + N_inv
            g0 += cdm.γ[j] * E0 * NiOverNj * oneOverSqrtEγE[j]
        end

        return g0 * -AE
    end
end



# Hole model parametrization from [1] equations (22)-(26), adjusted
@fastmath Λ(vrel::T) where {T} = T(0.75 * (1.0 - vrel))

@fastmath function Ω(vrel::T, ::Type{HPGe})::T where {T}
    p1::T = -0.29711
    p2::T = -1.12082
    p3::T = 3.83929
    p4::T = 4.80825
    x::T = 1 - vrel
    p1 * x + p2 * x^2 + p3 * x^3 + p4 * x^4
end

@fastmath function Ω(vrel::T, ::Type{Si})::T where {T}
    p1::T = -0.30565
    p2::T = -1.19650
    p3::T = 4.69001
    p4::T = -7.00635
    x::T = 1 - vrel
    p1 * x + p2 * x^2 + p3 * x^3 + p4 * x^4
end

@fastmath function getVh(fv::SVector{3,T}, cdm::ADLChargeDriftModel{T, M}, Emag_threshold::T = T(1e-5))::SVector{3,T} where {T <: SSDFloat, M}
    @inbounds begin
        Emag::T = norm(fv)
        Emag_inv::T = inv(Emag)

        if Emag < Emag_threshold return SVector{3,T}(0, 0, 0) end

        f::NTuple{4,T} = scale_to_given_temperature(Emag, cdm.temperaturemodel)
        f100h::T = f[3]
        f111h::T = f[4]

        V100h::T = Vl(Emag, cdm.holes.axis100) * f100h
        V111h::T = Vl(Emag, cdm.holes.axis111) * f111h
        
        tmp::SVector{3,T} = cdm.crystal_orientation * fv

        vrel::T = V111h / V100h
        Λvrel::T = Λ(vrel)
        Ωvrel::T = Ω(vrel, M)
        θ0::T = acos(tmp[3] / Emag)
        φ0::T = atan(tmp[2], tmp[1])
        sθ0::T, cθ0::T = sincos(θ0)
        s2θ0::T, c2θ0::T = sincos(2*θ0)
        s2φ0::T, c2φ0::T = sincos(2*φ0)

        vr::T = V100h * ( 1 - Λvrel * (sθ0^4 * s2φ0^2 + s2θ0^2) )
        vΩ::T = V100h * Ωvrel * 2 * (sθ0^3 * cθ0 * s2φ0^2 + s2θ0 * c2θ0)
        vφ::T = V100h * Ωvrel * 2 * sθ0^3 * s2φ0 * c2φ0

        Ry::RotY{T} = RotY{T}(θ0)
        Rz::SMatrix{3,3,T,9} = transpose(cdm.crystal_orientation) * RotZ{T}(φ0)
        Rz * (Ry * @SVector T[vΩ, vφ, vr])
    end
end
