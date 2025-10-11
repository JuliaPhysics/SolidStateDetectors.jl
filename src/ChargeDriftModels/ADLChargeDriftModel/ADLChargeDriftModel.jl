# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

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

# Longitudinal drift velocity formula
@fastmath function Vl(Emag::T, params::VelocityParameters{T})::T where {T <: SSDFloat}
    params.mu0 * Emag / (1 + (Emag / params.E0)^params.beta)^(1 / params.beta) - params.mun * Emag
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
include("ADL2006ChargeDriftModel.jl")
include("ADL2016ChargeDriftModel.jl")


# Electron model parametrization from Mihailescu (2000)
# https://doi.org/10.1016/S0168-9002(99)01286-3
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

getVe(fv::SVector{3, T}, cdm::ADLChargeDriftModel{T}, ::CartesianPoint{T}, Emag_threshold::T = T(1e-5)) where {T} = getVe(fv, cdm, Emag_threshold)
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

# Electron model parameterization from Bruyneel (2006), equations (1) - (8)
ν(E::T, params::ADLParameters{T}) where {T <: SSDFloat} = E ^ (params.Γ0 + params.Γ1 * log(E / params.Γ2))
ν(E::SVector{3,T}, params::ADLParameters{T}) where {T <: SSDFloat} = ν(norm(E), params)

@fastmath function SolidStateDetectors.getVe(fv::SVector{3,T}, cdm::ADL2016ChargeDriftModel{T,M}, Emag_threshold::T = T(1e-5))::SVector{3,T} where {T <: SSDFloat, M}
    @inbounds begin
        Γ0::T = sqrt((2 * cdm.parameters.mt_inv + cdm.parameters.ml_inv)/3)
        sqrtγ = sqrt.(cdm.γ)
        invν = broadcast(α -> inv(ν(α * fv, cdm.parameters)), sqrtγ)
        
        v::SVector{3,T} = SVector{3,T}([0,0,0])
        for j in eachindex(cdm.γ)
            Ei::SVector{3,T} = sqrtγ[j] * fv
            Ei_mag::T = norm(Ei)
            μ::T = Vl(Ei_mag / Γ0, cdm.electrons) / (Γ0 * Ei_mag)
            v -= invν[j] / sum(invν) * μ * cdm.γ[j] * fv
        end

        return v
    end
end


# Hole model parametrization from Bruyneel (2006) equations (22)-(26), adjusted
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

getVh(fv::SVector{3,T}, cdm::Union{ADLChargeDriftModel{T}, ADL2016ChargeDriftModel{T}}, ::CartesianPoint{T}, Emag_threshold::T = T(1e-5)) where {T <: SSDFloat} = getVh(fv, cdm, Emag_threshold)
@fastmath function getVh(fv::SVector{3,T}, cdm::Union{ADLChargeDriftModel{T, M}, ADL2016ChargeDriftModel{T, M}}, Emag_threshold::T = T(1e-5))::SVector{3,T} where {T <: SSDFloat, M}
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
