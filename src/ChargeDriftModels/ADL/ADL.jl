struct VelocityParameters{T <: SSDFloat}
    mu0::T
    beta::T
    E0::T
    mun::T
end

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
include("LinearModel.jl")
include("BoltzmannModel.jl")
include("PowerLawModel.jl")
include("VacuumModel.jl")

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


"""
    ADLChargeDriftModel{T <: SSDFloat, M <: AbstractDriftMaterial, N, TM <: AbstractTemperatureModel{T}} <: AbstractChargeDriftModel{T}

# Fields
- `electrons::CarrierParameters{T}`
- `holes::CarrierParameters{T}`
- `crystal_orientation::SMatrix{3,3,T,9}`
- `γ::SVector{N,SMatrix{3,3,T,9}}`
- `parameters::ADLParameters{T}`
- `temperaturemodel::TM`
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


ADLChargeDriftModel{T}(args...; kwargs...) where {T <: SSDFloat} = ADLChargeDriftModel(args..., T=T, kwargs...)
function ADLChargeDriftModel(configfilename::Union{Missing, AbstractString} = missing; T::Type=Float32, material::Type{<:AbstractDriftMaterial} = HPGe,
                             temperature::Union{Missing, Real}= missing, phi110::Union{Missing, Real} = missing)::ADLChargeDriftModel{T}

    if ismissing(configfilename) configfilename = joinpath(get_path_to_example_config_files(), "ADLChargeDriftModel/drift_velocity_config.json") end
    if !ismissing(temperature) temperature = T(temperature) end  #if you give the temperature it will be used, otherwise read from config file

    config = JSON.parsefile(configfilename)

    #mu0 in m^2 / ( V * s )
    #beta dimensionless
    #E0 in V / m
    #mun in m^2 / ( V * s )

    e100mu0::T  = config["drift"]["velocity"]["parameters"]["e100"]["mu0"]
    e100beta::T = config["drift"]["velocity"]["parameters"]["e100"]["beta"]
    e100E0::T   = config["drift"]["velocity"]["parameters"]["e100"]["E0"]
    e100mun::T  = config["drift"]["velocity"]["parameters"]["e100"]["mun"]
    e111mu0::T  = config["drift"]["velocity"]["parameters"]["e111"]["mu0"]
    e111beta::T = config["drift"]["velocity"]["parameters"]["e111"]["beta"]
    e111E0::T   = config["drift"]["velocity"]["parameters"]["e111"]["E0"]
    e111mun::T  = config["drift"]["velocity"]["parameters"]["e111"]["mun"]
    h100mu0::T  = config["drift"]["velocity"]["parameters"]["h100"]["mu0"]
    h100beta::T = config["drift"]["velocity"]["parameters"]["h100"]["beta"]
    h100E0::T   = config["drift"]["velocity"]["parameters"]["h100"]["E0"]
    h111mu0::T  = config["drift"]["velocity"]["parameters"]["h111"]["mu0"]
    h111beta::T = config["drift"]["velocity"]["parameters"]["h111"]["beta"]
    h111E0::T   = config["drift"]["velocity"]["parameters"]["h111"]["E0"]

    e100 = VelocityParameters{T}(e100mu0, e100beta, e100E0, e100mun)
    e111 = VelocityParameters{T}(e111mu0, e111beta, e111E0, e111mun)
    h100 = VelocityParameters{T}(h100mu0, h100beta, h100E0, 0)
    h111 = VelocityParameters{T}(h111mu0, h111beta, h111E0, 0)
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
            RotZ{T}(- π/4 - config["phi110"] )
        else
            transpose(parse_rotation_matrix(T, config, u"rad")) # replace u"rad" with the units in the config file
        end
    else
        RotZ{T}(- π/4 - phi110)
    end

    γ = setup_γj(crystal_orientation, parameters, material)

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


#= This function should never be called! ADLChargeDriftModel should be saved in the right format to the simulation!!
function getVe(fv::SVector{3, T}, cdm::ADLChargeDriftModel{<:Any,M,N,TM}, Emag_threshold::T = T(1e-5))::SVector{3, T} where {T <: SSDFloat, M, N, TM}
    @warn "ADLChargeDriftModel does not have the same precision type as the electric field vector."
    cdmT = ADLChargeDriftModel{T,M,N,TM}(cdm)
    getVe(fv, cdmT, Emag_threshold)
end
=#

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

        f::NTuple{4,T} = scale_to_given_temperature(cdm.temperaturemodel)
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

#= This function should never be called! ADLChargeDriftModel should be saved in the right format to the simulation!!
function getVh(fv::SVector{3,T}, cdm::ADLChargeDriftModel{<:Any,M,N,TM}, Emag_threshold::T = T(1e-5))::SVector{3,T} where {T <: SSDFloat, M, N, TM}
    @warn "ADLChargeDriftModel does not have the same precision type as the electric field vector."
    cdmT = ADLChargeDriftModel{T,M,N,TM}(cdm)
    getVh(fv, cdmT, Emag_threshold)
end
=#

@fastmath function getVh(fv::SVector{3,T}, cdm::ADLChargeDriftModel{T, M}, Emag_threshold::T = T(1e-5))::SVector{3,T} where {T <: SSDFloat, M}
    @inbounds begin
        Emag::T = norm(fv)
        Emag_inv::T = inv(Emag)

        if Emag < Emag_threshold return SVector{3,T}(0, 0, 0) end

        f::NTuple{4,T} = scale_to_given_temperature(cdm.temperaturemodel)
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


print(io::IO, tm::VacuumModel{T}) where {T <: SSDFloat} = print(io, "No temperature model defined")
println(io::IO, tm::VacuumModel) = print(io, tm)

print(io::IO, tm::BoltzmannModel{T}) where {T <: SSDFloat} = print(io, "BoltzmannModel{$T}")
function println(io::IO, tm::BoltzmannModel{T}) where {T <: SSDFloat}
    println("\n________BoltzmannModel________")
    println("Fit function: p1 + p2 exp(-p3/T)\n")
    println("---Temperature settings---")
    println("Crystal temperature:   \t $(tm.temperature)")
    println("Reference temperature: \t $(tm.reftemperature)\n")

    println("---Fitting parameters---")
    println("   \te100      \te111      \th100      \th111")
    println("p1 \t$(tm.p1e100)   \t$(tm.p1e111)   \t$(tm.p1h100)   \t$(tm.p1h111)")
    println("p2 \t$(tm.p2e100)   \t$(tm.p2e111)   \t$(tm.p2h100)   \t$(tm.p2h111)")
    println("p3 \t$(tm.p3e100)   \t$(tm.p3e111)   \t$(tm.p3h100)   \t$(tm.p3h111)")
end

print(io::IO, tm::LinearModel{T}) where {T <: SSDFloat} = print(io, "LinearModel{$T}")
function println(io::IO, tm::LinearModel{T}) where {T <: SSDFloat}
    println("\n________LinearModel________")
    println("Fit function: p1 + p2 * T\n")
    println("---Temperature settings---")
    println("Crystal temperature:  \t$(tm.temperature)")
    println("Reference temperature:\t$(tm.reftemperature)\n")

    println("---Fitting parameters---")
    println("   \te100      \te111      \th100      \th111")
    println("p1 \t$(tm.p1e100)   \t$(tm.p1e111)   \t$(tm.p1h100)   \t$(tm.p1h111)")
    println("p2 \t$(tm.p2e100)   \t$(tm.p2e111)   \t$(tm.p2h100)   \t$(tm.p2h111)")
end

print(io::IO, tm::PowerLawModel{T}) where {T <: SSDFloat} = print(io, "PowerLawModel{$T}")
function println(io::IO, tm::PowerLawModel{T}) where {T <: SSDFloat}
    println("\n________PowerLawModel________")
    println("Fit function: p1 * T^(3/2)\n")
    println("---Temperature settings---")
    println("Crystal temperature:   \t $(tm.temperature)")
    println("Reference temperature: \t $(tm.reftemperature)\n")

    println("---Fitting parameters---")
    println("   \te100      \te111      \th100      \th111")
    println("p1 \t$(tm.p1e100)   \t$(tm.p1e111)   \t$(tm.p1h100)   \t$(tm.p1h111)")
end

show(io::IO, tm::AbstractTemperatureModel) = print(io, tm) 
show(io::IO,::MIME"text/plain", tm::AbstractTemperatureModel) = show(io, tm)