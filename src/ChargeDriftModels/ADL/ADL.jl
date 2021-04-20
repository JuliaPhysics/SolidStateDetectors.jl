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

struct MassParameters{T <: SSDFloat}
    ml::T
    mt::T
end

# Temperature models
include("LinearModel.jl")
include("BoltzmannModel.jl")
include("PowerLawModel.jl")
include("VacuumModel.jl")

# Electron model parametrization from [3]
@fastmath function γj(j::Integer, phi110::T, γ0::SArray{Tuple{3,3},T,2,9}, ::Type{HPGe})::SArray{Tuple{3,3},T,2,9} where {T <: SSDFloat}
    tmp::T = 2 / 3
    a::T = acos(sqrt(tmp))
    Rx::SArray{Tuple{3,3},T,2,9} = SMatrix{3,3,T}(1, 0, 0, 0, cos(a), sin(a), 0, -sin(a), cos(a))
    b::T = phi110 + (j - 1) * T(π) / 2
    Rzj::SArray{Tuple{3,3},T,2,9} = SMatrix{3,3,T}(cos(b), -sin(b), 0, sin(b), cos(b), 0, 0, 0, 1)
    Rj = Rx * Rzj
    transpose(Rj) * γ0 * Rj
end

@fastmath function γj(j::Integer, phi110::T, γ0::SArray{Tuple{3,3},T,2,9}, ::Type{Si})::SArray{Tuple{3,3},T,2,9} where {T <: SSDFloat}
    b::T = phi110 + T(π/4) * (-1)^j
    Rj::SArray{Tuple{3,3},T,2,9} = (j == 3 ? SMatrix{3,3,T}(1, 0, 0, 0, 0, 1, 0, -1, 0)  : I) * SMatrix{3,3,T}(cos(b), -sin(b), 0, sin(b), cos(b), 0, 0, 0, 1) 
    transpose(Rj) * γ0 * Rj
end

@fastmath function setup_gamma_matrices(phi110::T, masses::MassParameters{T}, material::Type{HPGe})::SVector{4, SMatrix{3,3,T}} where {T <: SSDFloat}
    ml::T = masses.ml
    mt::T = masses.mt
    γ0 = SMatrix{3,3,T}(1. / mt, 0, 0, 0, 1. / ml, 0, 0, 0, 1. / mt)
    SVector{4, SArray{Tuple{3,3},T,2,9}}([γj(i, phi110, γ0, material) for i in Base.OneTo(4)]...)
end

@fastmath function setup_gamma_matrices(phi110::T, masses::MassParameters{T}, material::Type{Si})::SVector{3, SMatrix{3,3,T}} where {T <: SSDFloat}
    ml::T = masses.ml
    mt::T = masses.mt
    γ0 = SMatrix{3,3,T}(1. / mt, 0, 0, 0, 1. / ml, 0, 0, 0, 1. / mt)
    SVector{3, SArray{Tuple{3,3},T,2,9}}([γj(i, phi110, γ0, material) for i in Base.OneTo(3)]...)
end

# Longitudinal drift velocity formula
@fastmath function Vl(Emag::T, params::VelocityParameters{T})::T where {T <: SSDFloat}
    params.mu0 * Emag / (1 + (Emag / params.E0)^params.beta)^(1 / params.beta) - params.mun * Emag
end


"""
    ADLChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel

# Fields
- `material::Type{<:AbstractDriftMaterial}`
- `electrons::CarrierParameters{T}`
- `holes::CarrierParameters{T}`
- `masses::MassParameters{T}`
- `phi110::T
- `gammas::Vector{SMatrix{3,3,T}}`
- `temperaturemodel::AbstractTemperatureModel{T}`
"""
struct ADLChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
    material::Type{<:AbstractDriftMaterial}
    electrons::CarrierParameters{T}
    holes::CarrierParameters{T}
    masses::MassParameters{T}
    phi110::T
    gammas::Vector{SMatrix{3,3,T}}
    temperaturemodel::AbstractTemperatureModel{T}
end

function ADLChargeDriftModel{T}(chargedriftmodel::ADLChargeDriftModel)::ADLChargeDriftModel{T} where {T <: SSDFloat}
    cdmf64 = chargedriftmodel
    e100 = VelocityParameters{T}(cdmf64.electrons.axis100.mu0, cdmf64.electrons.axis100.beta, cdmf64.electrons.axis100.E0, cdmf64.electrons.axis100.mun)
    e111 = VelocityParameters{T}(cdmf64.electrons.axis111.mu0, cdmf64.electrons.axis111.beta, cdmf64.electrons.axis111.E0, cdmf64.electrons.axis111.mun)
    h100 = VelocityParameters{T}(cdmf64.holes.axis100.mu0, cdmf64.holes.axis100.beta, cdmf64.holes.axis100.E0, cdmf64.holes.axis100.mun)
    h111 = VelocityParameters{T}(cdmf64.holes.axis111.mu0, cdmf64.holes.axis111.beta, cdmf64.holes.axis111.E0, cdmf64.holes.axis111.mun)
    electrons = CarrierParameters{T}(e100, e111)
    holes     = CarrierParameters{T}(h100, h111)
    ml::T     = cdmf64.masses.ml
    mt::T     = cdmf64.masses.mt
    masses    = MassParameters{T}(ml, mt)
    phi110::T = cdmf64.phi110
    gammas = Vector{SArray{Tuple{3,3},T,2,9}}( cdmf64.gammas )
    temperaturemodel::AbstractTemperatureModel{T} = cdmf64.temperaturemodel
    material = cdmf64.material
    ADLChargeDriftModel{T}(material, electrons, holes, masses, phi110, gammas, temperaturemodel)
end

function ADLChargeDriftModel(configfilename::Union{Missing, AbstractString} = missing; T::Type=Float32, material::Type{<:AbstractDriftMaterial} = HPGe,
                             temperature::Union{Missing, Real}= missing, phi110::Union{Missing, Real} = missing)::ADLChargeDriftModel{T}

    if ismissing(configfilename) configfilename = joinpath(@__DIR__, "drift_velocity_config.json") end
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
        (config["material"] == "Si") && (material = Si)
        (config["material"] == "HPGe") && (material = HPGe)
    end

    masses = if "masses" in keys(config["drift"])
        ml = T(config["drift"]["masses"]["ml"])
        mt = T(config["drift"]["masses"]["mt"])
        MassParameters{T}(ml, mt)
    else
        ml = T(material_properties[Symbol(material)].ml)
        mt = T(material_properties[Symbol(material)].mt)
        MassParameters{T}(ml, mt)
    end

    ismissing(phi110) ? phi110 = T(config["phi110"]) : phi110 = T(phi110)  #if you give the angle of the 110 axis it will be used, otherwise read from config file

    gammas = setup_gamma_matrices(phi110, masses, material)

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

    return ADLChargeDriftModel{T}(material, electrons, holes, masses, phi110, gammas, temperaturemodel)
end



# This should never be called! ADLChargeDriftModel should be saved in the right format to the simulation!!
function getVe(fv::SVector{3, T}, cdm::ADLChargeDriftModel, Emag_threshold::T = T(1e-5))::SVector{3, T} where {T <: SSDFloat}
    @warn "ADLChargeDriftModel does not have the same precision type as the electric field vector."
    cdmT = ADLChargeDriftModel{T}(cdm)
    getVe(fv, cdmT, Emag_threshold)
end

@inline get_AE_and_RE(cdm::ADLChargeDriftModel{T}, V100e::T, V111e::T) where {T} = _get_AE_and_RE(cdm, cdm.material, V100e, V111e)

function _get_AE_and_RE(cdm::ADLChargeDriftModel{T}, ::Type{HPGe}, V100e::T, V111e::T)::Tuple{T,T} where{T <: SSDFloat}
    mli::T = 1. / cdm.masses.ml
    mti::T = 1. / cdm.masses.mt
    B::T = 1/3 * sqrt(8*mti + mli)

    Gamma0::T = sqrt(2/3*mti + 1/3*mli)
    Gamma1::T = (-4*sqrt(mli) -4/3 * B)/(sqrt(mli) - B)^2
    Gamma2::T = Gamma1*(sqrt(mli) + 3*B)/(-4)

    AE::T = V100e / Gamma0
    RE::T = Gamma1 * V111e / AE + Gamma2
    
    AE, RE
end

function _get_AE_and_RE(cdm::ADLChargeDriftModel{T}, ::Type{Si}, V100e::T, V111e::T)::Tuple{T,T} where{T <: SSDFloat}
    
    mli::T = 1. / cdm.masses.ml
    mti::T = 1. / cdm.masses.mt

    Gamma0::T = sqrt(2/3*mti + 1/3*mli)
    Gamma1::T = -1.5 * (sqrt(mti) + 2*sqrt(mli)) / (sqrt(mli) - sqrt(mti))^2
    Gamma2::T = 0.5 * (sqrt(mli) + 2*sqrt(mti)) * (2*sqrt(mli) + sqrt(mti)) / (sqrt(mli) - sqrt(mti))^2
    
    AE::T = V111e / Gamma0
    RE::T = Gamma1 * V100e / AE + Gamma2
    
    AE, RE
end


@fastmath function getVe(fv::SVector{3, T}, cdm::ADLChargeDriftModel{T}, Emag_threshold::T = T(1e-5))::SVector{3, T} where {T <: SSDFloat}
    @inbounds begin
        Emag::T = norm(fv)
        Emag_inv::T = inv(Emag)

        if Emag < Emag_threshold return SVector{3,T}(0, 0, 0) end

        f::NTuple{4,T} = scale_to_given_temperature(cdm.temperaturemodel)
        f100e::T = f[1]
        f111e::T = f[2]
        V100e::T = Vl(Emag, cdm.electrons.axis100) * f100e
        V111e::T = Vl(Emag, cdm.electrons.axis111) * f111e
        
        AE::T, RE::T = get_AE_and_RE(cdm, V100e, V111e)

        e0 = SVector{3, T}(fv * Emag_inv)
        
        N::Int = length(cdm.gammas) # number of degenerate conduction band minima
        oneOverSqrtEgE::SVector{N, T} = [1 / sqrt( cdm.gammas[j] * e0 ⋅ e0 ) for j in eachindex(cdm.gammas)] # setup.gammas[j] * e0 ⋅ e0 -> causes allocations

        g0 = MVector{3, T}(0,0,0)
        for j in eachindex(cdm.gammas)
            NiOverNj::T = RE * (oneOverSqrtEgE[j] / sum(oneOverSqrtEgE) - T(1/N)) + T(1/N)
            g0 += cdm.gammas[j] * e0 * NiOverNj * oneOverSqrtEgE[j]
        end

        return g0 * -AE
    end
end



# Hole model parametrization from [1] equations (22)-(26)
@fastmath function k0func(vrel::T, ::Type{HPGe})::T where {T <: SSDFloat}
    p0::T = 9.2652
    p1::T = 26.3467
    p2::T = 29.6137
    p3::T = 12.3689
    return @fastmath p0 - p1 * vrel + p2 * vrel^2 - p3 * vrel^3
end

@fastmath function k0func(vrel::T, ::Type{Si})::T where {T <: SSDFloat}
    p0::T = 14.2579
    p1::T = 42.3839
    p2::T = 47.4626
    p3::T = 19.1738
    return @fastmath p0 - p1 * vrel + p2 * vrel^2 - p3 * vrel^3
end

@fastmath function lambda(k0::T, ::Type{HPGe})::T where {T <: SSDFloat}
    p0::T = -0.01322
    p1::T = 0.41145
    p2::T = 0.23657
    p3::T = 0.04077
    return @fastmath p0 * k0 + p1 * k0^2 - p2 * k0^3 + p3 * k0^4
end

@fastmath function lambda(k0::T, ::Type{Si})::T where {T <: SSDFloat}
    p0::T = -0.01687
    p1::T = 0.34789
    p2::T = 0.19211
    p3::T = 0.03171
    return @fastmath p0 * k0 + p1 * k0^2 - p2 * k0^3 + p3 * k0^4
end

@fastmath function omega(k0::T, ::Type{HPGe})::T where {T <: SSDFloat}
    p0::T = 0.006550
    p1::T = 0.19946
    p2::T = 0.09859
    p3::T = 0.01559
    return @fastmath p0 * k0 - p1 * k0^2 + p2 * k0^3 - p3 * k0^4
end

@fastmath function omega(k0::T, ::Type{Si})::T where {T <: SSDFloat}
    p0::T = 0.00670
    p1::T = 0.16390
    p2::T = 0.07575
    p3::T = 0.01091
    return @fastmath p0 * k0 - p1 * k0^2 + p2 * k0^3 - p3 * k0^4
end


# This function should never be called!
function getVh(fv::SVector{3,T}, cdm::ADLChargeDriftModel, Emag_threshold::T = T(1e-5))::SVector{3,T} where {T <: SSDFloat}
    @warn "ADLChargeDriftModel does not have the same precision type as the electric field vector."
    cdmT = ADLChargeDriftModel{T}(cdm)
    getVh(fv, cdmT, Emag_threshold)
end

@inline function getVh(fv::SVector{3,T}, cdm::ADLChargeDriftModel{T}, Emag_threshold::T = T(1e-5))::SVector{3,T} where {T <: SSDFloat}
    getVh(fv, cdm, cdm.material, Emag_threshold)
end

@fastmath function getVh(fv::SVector{3,T}, cdm::ADLChargeDriftModel{T}, material::Type{<:AbstractDriftMaterial}, Emag_threshold::T = T(1e-5))::SVector{3,T} where {T <: SSDFloat}
    @inbounds begin
        Emag::T = norm(fv)
        Emag_inv::T = inv(Emag)

        if Emag < Emag_threshold return SVector{3,T}(0, 0, 0) end

        f::NTuple{4,T} = scale_to_given_temperature(cdm.temperaturemodel)
        f100h::T = f[3]
        f111h::T = f[4]

        V100h::T = Vl(Emag, cdm.holes.axis100) * f100h
        V111h::T = Vl(Emag, cdm.holes.axis111) * f111h

        b::T = π / 4 + cdm.phi110
        Rz = SMatrix{3, 3, T}(cos(b), -sin(b), 0, sin(b), cos(b), 0, 0, 0, 1)
        a = Rz * fv

        k0::T = k0func(V111h / V100h, material)
        Λk0::T = lambda(k0, material)  # alternative without k0: 0.75 * (V100h - V11h) / V100h
        Ωk0::T = omega(k0, material)
        θ0::T = acos(a[3] / Emag)
        φ0::T = atan(a[2], a[1])

        vtmp = MVector{3, T}(0, 0, 0)
        vtmp[3] = V100h * ( 1 - Λk0 * (sin(θ0)^4 * sin(2 * φ0)^2 + sin(2 * θ0)^2) )
        vtmp[1] = V100h * Ωk0 * (2 * sin(θ0)^3 * cos(θ0) * sin(2 * φ0)^2 + sin(4 * θ0))
        vtmp[2] = V100h * Ωk0 * sin(θ0)^3 * sin(4 * φ0)

        Ry = SMatrix{3, 3, T}(cos(θ0),0,-sin(θ0), 0,1,0, sin(θ0),0,cos(θ0))
        b = φ0 + π / 4 + cdm.phi110
        Rz = SMatrix{3, 3, T}(cos(b), sin(b), 0, -sin(b), cos(b), 0, 0, 0, 1)
        Rz * (Ry * vtmp)
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