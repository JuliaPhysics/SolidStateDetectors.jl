abstract type AbstractChargeDriftModels end


"""
    VacuumChargeDriftModel <: AbstractChargeDriftModels
"""
struct VacuumChargeDriftModel <: AbstractChargeDriftModels end

function get_electron_drift_field(ef::Array{SVector{3,T},3}, ::VacuumChargeDriftModel)::Array{SVector{3,T},3} where {T<:AbstractFloat}
    return -ef
end
function get_hole_drift_field(ef::Array{SVector{3,T},3}, ::VacuumChargeDriftModel)::Array{SVector{3,T},3} where {T<:AbstractFloat}
    return ef
end



#################################
### Start: ADL Charge Drift Model

struct VelocityParameters{T}
    mu0::T
    beta::T
    E0::T
    mun::T
end

struct CarrierParameters{T}
    axis100::VelocityParameters{T}
    axis111::VelocityParameters{T}
end

# Electron model parametrization from [3]
function γj(j::Integer, γ0::SArray{Tuple{3,3},T,2,9})::SArray{Tuple{3,3},T,2,9} where {T <: AbstractFloat}
    tmp::T = 2 / 3
    a::T = acos(sqrt(tmp))
    Rx::SArray{Tuple{3,3},T,2,9} = SMatrix{3,3,T}(1, 0, 0, 0, cos(a), sin(a), 0, -sin(a), cos(a))
    b::T = (j - 1) * T(π) / 2
    Rzj::SArray{Tuple{3,3},T,2,9} = SMatrix{3,3,T}(cos(b), sin(b), 0, -sin(b), cos(b), 0, 0, 0, 1)
    Rj = Rx * Rzj
    transpose(Rj) * γ0 * Rj
end

function setup_gamma_matrices(phi110::T)::SVector{4, SMatrix{3,3,T}} where {T <: AbstractFloat}
    ml::T = 1.64
    mt::T = 0.0819
    γ0 = SMatrix{3,3,T}(1 / mt, 0, 0, 0, 1 / ml, 0, 0, 0, 1 / mt)
    SVector{4, SArray{Tuple{3,3},T,2,9}}(γj(1, γ0), γj(2, γ0), γj(3, γ0), γj(4, γ0))
end

# Longitudinal drift velocity formula
@fastmath function Vl(Emag::T, params::VelocityParameters{T})::T where {T <: AbstractFloat}
    params.mu0 * Emag / (1 + (Emag / params.E0)^params.beta)^(1 / params.beta) - params.mun * Emag
end


"""
    ADLChargeDriftModel{T <: AbstractFloat} <: AbstractChargeDriftModels

# Fields
- `electrons::CarrierParameters{T}`
- `holes::CarrierParameters{T}`
- `phi110::T
- `gammas::SVector{4, SMatrix{3,3,T}}`
"""
struct ADLChargeDriftModel{T <: AbstractFloat} <: AbstractChargeDriftModels 
    electrons::CarrierParameters{T}
    holes::CarrierParameters{T}
    phi110::T
    gammas::SVector{4, SMatrix{3,3,T}}
end

function ADLChargeDriftModel(configfilename::Union{Missing, AbstractString} = missing; T::Type=Float64)::ADLChargeDriftModel{T} 
    if ismissing(configfilename) configfilename = joinpath(@__DIR__, "drift_velocity_config.json") end
    
    config = JSON.parsefile(configfilename)

    #mu0 in m^2 / ( V * s )
    e100mu0::T  = config["drift"]["velocity"]["parameters"]["e100"]["mu0"]
    #beta dimensionless
    e100beta::T = config["drift"]["velocity"]["parameters"]["e100"]["beta"]
    #E0 in V / m
    e100E0::T   = config["drift"]["velocity"]["parameters"]["e100"]["E0"]
    #mun in m^2 / ( V * s )
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

    phi110::T = config["phi110"]

    gammas = setup_gamma_matrices(phi110)

    return ADLChargeDriftModel{T}(electrons, holes, phi110, gammas)
end

function get_electron_drift_field(ef::Array{SVector{3,T},3}, chargedriftmodel::ADLChargeDriftModel)::Array{SVector{3,T},3} where {T<:AbstractFloat}
    df = Array{SVector{3,T}, 3}(undef, size(ef))

    cdm = begin
        cdmf64 = chargedriftmodel
        e100 = VelocityParameters{T}(cdmf64.electrons.axis100.mu0, cdmf64.electrons.axis100.beta, cdmf64.electrons.axis100.E0, cdmf64.electrons.axis100.mun)
        e111 = VelocityParameters{T}(cdmf64.electrons.axis111.mu0, cdmf64.electrons.axis111.beta, cdmf64.electrons.axis111.E0, cdmf64.electrons.axis111.mun)
        h100 = VelocityParameters{T}(cdmf64.holes.axis100.mu0, cdmf64.holes.axis100.beta, cdmf64.holes.axis100.E0, cdmf64.holes.axis100.mun)
        h111 = VelocityParameters{T}(cdmf64.holes.axis111.mu0, cdmf64.holes.axis111.beta, cdmf64.holes.axis111.E0, cdmf64.holes.axis111.mun)
        electrons = CarrierParameters{T}(e100, e111)
        holes     = CarrierParameters{T}(h100, h111)  
        phi110::T = cdmf64.phi110
        gammas = SVector{4, SArray{Tuple{3,3},T,2,9}}( cdmf64.gammas )
        ADLChargeDriftModel{T}(electrons, holes, phi110, gammas)
    end

    function getVe(fv::SVector{3,T}, gammas::SVector{4, SMatrix{3,3,T}})::SVector{3,T} where {T <: AbstractFloat}
        @inbounds begin
            Emag::T = norm(fv)
            Emag_inv::T = inv(Emag)

            Emag_threshold::T = 1e-5

            tmp0::T = 2.888470213
            tmp1::T = -1.182108256
            tmp2::T = 3.160660533
            tmp3::T = 0.25

            if Emag < Emag_threshold
                return SVector{3,T}(0, 0, 0)
            end

            V100e::T = Vl(Emag, cdm.electrons.axis100)
            V111e::T = Vl(Emag, cdm.electrons.axis111)

            AE::T = V100e / tmp0
            RE::T = tmp1 * V111e / AE + tmp2

            e0 = SVector{3, T}(fv * Emag_inv)

            oneOverSqrtEgE::SVector{4, T} = [1 / sqrt( gammas[j] * e0 ⋅ e0 ) for j in eachindex(1:4)] # setup.gammas[j] * e0 ⋅ e0 -> causes allocations

            g0 = MVector{3, T}(0,0,0)
            for j in eachindex(1:4)
                NiOverNj::T = RE * (oneOverSqrtEgE[j] / sum(oneOverSqrtEgE) - tmp3) + tmp3
                g0 += gammas[j] * e0 * NiOverNj * oneOverSqrtEgE[j]
            end

            return g0 * -AE
        end
    end
    for i in eachindex(df)
        @inbounds df[i] = getVe(ef[i], cdm.gammas)
    end

    return df
end

# Hole model parametrization from [1] equations (22)-(26)
function k0func(vrel::T)::T where {T <: AbstractFloat}
    p0::T = 9.2652
    p1::T = 26.3467
    p2::T = 29.6137
    p3::T = 12.3689
    return @fastmath p0 - p1 * vrel + p2 * vrel^2 - p3 * vrel^3
end

@fastmath function lambda(k0::T)::T where {T <: AbstractFloat}
    p0::T = -0.01322
    p1::T = 0.41145
    p2::T = 0.23657
    p3::T = 0.04077
    return @fastmath p0 * k0 + p1 * k0^2 - p2 * k0^3 + p3 * k0^4
end

@fastmath function omega(k0::T)::T where {T <: AbstractFloat}
    p0::T = 0.006550
    p1::T = 0.19946
    p2::T = 0.09859
    p3::T = 0.01559
    return @fastmath p0 * k0 - p1 * k0^2 + p2 * k0^3 - p3 * k0^4
end


function get_hole_drift_field(ef::Array{SVector{3,T},3}, chargedriftmodel::ADLChargeDriftModel)::Array{SVector{3,T},3} where {T<:AbstractFloat}
    df = Array{SVector{3,T}, 3}(undef, size(ef))

    cdm = begin
        cdmf64 = chargedriftmodel
        e100 = VelocityParameters{T}(cdmf64.electrons.axis100.mu0, cdmf64.electrons.axis100.beta, cdmf64.electrons.axis100.E0, cdmf64.electrons.axis100.mun)
        e111 = VelocityParameters{T}(cdmf64.electrons.axis111.mu0, cdmf64.electrons.axis111.beta, cdmf64.electrons.axis111.E0, cdmf64.electrons.axis111.mun)
        h100 = VelocityParameters{T}(cdmf64.holes.axis100.mu0, cdmf64.holes.axis100.beta, cdmf64.holes.axis100.E0, cdmf64.holes.axis100.mun)
        h111 = VelocityParameters{T}(cdmf64.holes.axis111.mu0, cdmf64.holes.axis111.beta, cdmf64.holes.axis111.E0, cdmf64.holes.axis111.mun)
        electrons = CarrierParameters{T}(e100, e111)
        holes     = CarrierParameters{T}(h100, h111)  
        phi110::T = cdmf64.phi110
        gammas = SVector{4, SArray{Tuple{3,3},T,2,9}}( cdmf64.gammas )
        ADLChargeDriftModel{T}(electrons, holes, phi110, gammas)
    end

    Emag_threshold::T = 1e-5

    function getVh(fv::SVector{3,T}, Emag_threshold::T)::SVector{3,T} where {T <: AbstractFloat}
        @inbounds begin
            Emag::T = norm(fv)
            Emag_inv::T = inv(Emag)

            if Emag < Emag_threshold
                return SVector{3,T}(0, 0, 0)
            end
            V100h::T = Vl(Emag, cdm.holes.axis100)
            V111h::T = Vl(Emag, cdm.holes.axis111)

            b::T = -π / 4 - cdm.phi110
            Rz = SMatrix{3, 3, T}(cos(b), sin(b), 0, -sin(b), cos(b), 0, 0, 0, 1)
            a = Rz * fv

            theta0::T = acos(a[3] / Emag)
            phi0::T = atan(a[2], a[1])

            k0::T = k0func(V111h / V100h)

            vtmp = MVector{3, T}(0, 0, 0)
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

    for i in eachindex(df)
        @inbounds df[i] = getVh(ef[i], Emag_threshold)
    end

    return df
end


### END: ADL Charge Drift Model
###############################



