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
@fastmath function γj(j::Integer, phi110::T, γ0::SArray{Tuple{3,3},T,2,9})::SArray{Tuple{3,3},T,2,9} where {T <: SSDFloat}
    tmp::T = 2 / 3
    a::T = acos(sqrt(tmp))
    Rx::SArray{Tuple{3,3},T,2,9} = SMatrix{3,3,T}(1, 0, 0, 0, cos(a), sin(a), 0, -sin(a), cos(a))
    b::T = -phi110 + (j - 1) * T(π) / 2
    Rzj::SArray{Tuple{3,3},T,2,9} = SMatrix{3,3,T}(cos(b), sin(b), 0, -sin(b), cos(b), 0, 0, 0, 1)
    Rj = Rx * Rzj
    transpose(Rj) * γ0 * Rj
end

@fastmath function setup_gamma_matrices(phi110::T, masses::MassParameters{T})::SVector{4, SMatrix{3,3,T}} where {T <: SSDFloat}
    ml::T = masses.ml
    mt::T = masses.mt
    γ0 = SMatrix{3,3,T}(1. / mt, 0, 0, 0, 1. / ml, 0, 0, 0, 1. / mt)
    SVector{4, SArray{Tuple{3,3},T,2,9}}(γj(1, phi110, γ0), γj(2, phi110, γ0), γj(3, phi110, γ0), γj(4, phi110, γ0))
end

# Longitudinal drift velocity formula
@fastmath function Vl(Emag::T, params::VelocityParameters{T})::T where {T <: SSDFloat}
    params.mu0 * Emag / (1 + (Emag / params.E0)^params.beta)^(1 / params.beta) - params.mun * Emag
end


"""
    ADLChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel

# Fields
- `electrons::CarrierParameters{T}`
- `holes::CarrierParameters{T}`
- `masses::MassParameters{T}`
- `phi110::T
- `gammas::SVector{4, SMatrix{3,3,T}}`
- `temperaturemodel::AbstractTemperatureModel{T}`
"""
struct ADLChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
    electrons::CarrierParameters{T}
    holes::CarrierParameters{T}
    masses::MassParameters{T}
    phi110::T
    gammas::SVector{4, SMatrix{3,3,T}}
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
    gammas = SVector{4, SArray{Tuple{3,3},T,2,9}}( cdmf64.gammas )
    temperaturemodel::AbstractTemperatureModel{T} = cdmf64.temperaturemodel
    ADLChargeDriftModel{T}(electrons, holes, masses, phi110, gammas, temperaturemodel)
end

function ADLChargeDriftModel(configfilename::Union{Missing, AbstractString} = missing; T::Type=Float32,
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

    if "masses" in keys(config["drift"])
        ml = T(config["drift"]["masses"]["ml"])
        mt = T(config["drift"]["masses"]["mt"])
        masses = MassParameters{T}(ml, mt)
    else
        masses = MassParameters{T}(T(1.64),T(0.0819))
    end

    ismissing(phi110) ? phi110 = T(config["phi110"]) : phi110 = T(phi110)  #if you give the angle of the 110 axis it will be used, otherwise read from config file

    gammas = setup_gamma_matrices(phi110, masses)

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

    return ADLChargeDriftModel{T}(electrons, holes, masses, phi110, gammas, temperaturemodel)
end




function getVe(fv::SVector{3, T}, cdm::ADLChargeDriftModel, Emag_threshold::T = T(1e-5))::SVector{3, T} where {T <: SSDFloat}
    cdmT = ADLChargeDriftModel{T}(cdm)
    getVe(fv, cdmT, Emag_threshold)
end

@fastmath function getVe(fv::SVector{3, T}, cdm::ADLChargeDriftModel{T}, Emag_threshold::T = T(1e-5))::SVector{3, T} where {T <: SSDFloat}
    @inbounds begin
        Emag::T = norm(fv)
        Emag_inv::T = inv(Emag)

        if Emag < Emag_threshold
            return SVector{3,T}(0, 0, 0)
        end

        mli::T = 1. / cdm.masses.ml
        mti::T = 1. / cdm.masses.mt
        B::T = 1/3 * sqrt(8*mti + mli)

        Gamma0::T = sqrt(2/3*mti + 1/3*mli)
        Gamma1::T = (-4*sqrt(mli) -4/3 * B)/(sqrt(mli) - B)^2
        Gamma2::T = Gamma1*(sqrt(mli) + 3*B)/(-4)

        f::NTuple{4,T} = scale_to_given_temperature(cdm.temperaturemodel)
        f100e::T = f[1]
        f111e::T = f[2]

        V100e::T = Vl(Emag, cdm.electrons.axis100) * f100e
        V111e::T = Vl(Emag, cdm.electrons.axis111) * f111e

        AE::T = V100e / Gamma0
        RE::T = Gamma1 * V111e / AE + Gamma2

        e0 = SVector{3, T}(fv * Emag_inv)

        oneOverSqrtEgE::SVector{4, T} = [1 / sqrt( cdm.gammas[j] * e0 ⋅ e0 ) for j in eachindex(1:4)] # setup.gammas[j] * e0 ⋅ e0 -> causes allocations

        g0 = MVector{3, T}(0,0,0)
        for j in eachindex(1:4)
            NiOverNj::T = RE * (oneOverSqrtEgE[j] / sum(oneOverSqrtEgE) - T(0.25)) + T(0.25)
            g0 += cdm.gammas[j] * e0 * NiOverNj * oneOverSqrtEgE[j]
        end

        return g0 * -AE
    end
end



# Hole model parametrization from [1] equations (22)-(26)
@fastmath function k0func(vrel::T)::T where {T <: SSDFloat}
    p0::T = 9.2652
    p1::T = 26.3467
    p2::T = 29.6137
    p3::T = 12.3689
    return @fastmath p0 - p1 * vrel + p2 * vrel^2 - p3 * vrel^3
end

@fastmath function lambda(k0::T)::T where {T <: SSDFloat}
    p0::T = -0.01322
    p1::T = 0.41145
    p2::T = 0.23657
    p3::T = 0.04077
    return @fastmath p0 * k0 + p1 * k0^2 - p2 * k0^3 + p3 * k0^4
end

@fastmath function omega(k0::T)::T where {T <: SSDFloat}
    p0::T = 0.006550
    p1::T = 0.19946
    p2::T = 0.09859
    p3::T = 0.01559
    return @fastmath p0 * k0 - p1 * k0^2 + p2 * k0^3 - p3 * k0^4
end



function getVh(fv::SVector{3,T}, cdm::ADLChargeDriftModel, Emag_threshold::T = T(1e-5))::SVector{3,T} where {T <: SSDFloat}
    cdmT = ADLChargeDriftModel{T}(cdm)
    getVh(fv, cdmT, Emag_threshold)
end

@fastmath function getVh(fv::SVector{3,T}, cdm::ADLChargeDriftModel{T}, Emag_threshold::T = T(1e-5))::SVector{3,T} where {T <: SSDFloat}
    @inbounds begin
        Emag::T = norm(fv)
        Emag_inv::T = inv(Emag)

        if Emag < Emag_threshold
            return SVector{3,T}(0, 0, 0)
        end

        f::NTuple{4,T} = scale_to_given_temperature(cdm.temperaturemodel)
        f100h::T = f[3]
        f111h::T = f[4]

        V100h::T = Vl(Emag, cdm.holes.axis100) * f100h
        V111h::T = Vl(Emag, cdm.holes.axis111) * f111h

        b::T = -π / 4 - cdm.phi110
        Rz = SMatrix{3, 3, T}(cos(b), sin(b), 0, -sin(b), cos(b), 0, 0, 0, 1)
        a = Rz * fv

        theta0::T = acos(a[3] / Emag)
        phi0::T = atan(a[2], a[1])

        k0::T = k0func(V111h / V100h)

        vtmp = MVector{3, T}(0, 0, 0) ## from CITATION; The implementation here is correct, mistake in the CITATION
        vtmp[3] = V100h * ( 1 - lambda(k0) * (sin(theta0)^4 * sin(2 * phi0)^2 + sin(2 * theta0)^2) )
        vtmp[1] = V100h * omega(k0) * (2 * sin(theta0)^3 * cos(theta0) * sin(2 * phi0)^2 + sin(4 * theta0))
        vtmp[2] = V100h * omega(k0) * sin(theta0)^3 * sin(4 * phi0)

        Ry = SMatrix{3, 3, T}(cos(theta0),0,-sin(theta0), 0,1,0, sin(theta0),0,cos(theta0))
        b = phi0 + pi/4 + cdm.phi110
        Rz = SMatrix{3,3, T}(cos(b),sin(b),0, -sin(b),cos(b),0, 0,0,1)
        vtmp = Rz * (Ry * vtmp)

        return vtmp
    end
end




function print(io::IO, tm::VacuumModel{T}) where {T <: SSDFloat}
    print(io, "No temperature model defined")
end
println(io::IO, tm::VacuumModel) = print(io, tm)

function print(io::IO, tm::BoltzmannModel{T}) where {T <: SSDFloat}
    print(io, "BoltzmannModel{$T}")
end
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

function print(io::IO, tm::LinearModel{T}) where {T <: SSDFloat}
    print(io, "LinearModel{$T}")
end
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

function print(io::IO, tm::PowerLawModel{T}) where {T <: SSDFloat}
    print(io, "PowerLawModel{$T}")
end
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

