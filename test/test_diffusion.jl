using SolidStateDetectors
using SolidStateDetectors: Electron, Hole, AbstractChargeDriftModel
using LinearAlgebra
using Test

T = Float32

@testset "Test position-dependent diffusion" begin
    
    N = 10
    done = fill(false, N)
    Δt = T(1e-9)
    zs = 1:10
    current_pos = CartesianPoint{T}.(0, 0, zs)
    crystal_temperature = T(78)
    
    step_vectors = zeros(CartesianVector{T}, N)

    # Charge Drift model with zero mobility
    cdm0 = IsotropicChargeDriftModel{T}(0,0)
    SolidStateDetectors._add_fieldvector_diffusion!(step_vectors, done, Δt, cdm0, current_pos, crystal_temperature, Hole)
    @test all(iszero.(step_vectors))
    SolidStateDetectors._add_fieldvector_diffusion!(step_vectors, done, Δt, cdm0, current_pos, crystal_temperature, Electron)
    @test all(iszero.(step_vectors))
    
    # if the mobility is set to some position-independent value, then we know what the length of the expected step_vectors should be
    μ0 = T(5.0 + rand() * 5.0)
    cdmµ = IsotropicChargeDriftModel{T}(µ0, µ0)
    SolidStateDetectors._add_fieldvector_diffusion!(step_vectors, done, Δt, cdmµ, current_pos, crystal_temperature, Hole)
    D0 = μ0 * crystal_temperature * SolidStateDetectors.kB / SolidStateDetectors.elementary_charge
    @test all(isapprox.(norm.(step_vectors), sqrt(6 * D0 * Δt), rtol = 0.01))
    step_vectors = zeros(CartesianVector{T}, N) # reset step vectors
    
    # define a new charge drift model with position-dependent mobility
    struct PositionDependentChargeDriftModel{T} <: AbstractChargeDriftModel{T} end
    SolidStateDetectors.calculate_mobility(::PositionDependentChargeDriftModel{T}, pt::CartesianPoint{T}, ::Type{Hole}) where {T} = pt.z
    SolidStateDetectors._add_fieldvector_diffusion!(step_vectors, done, Δt, PositionDependentChargeDriftModel{T}(), current_pos, crystal_temperature, Hole)
    @test all(isapprox.(norm.(step_vectors).^2 ./ (6 * Δt) * SolidStateDetectors.elementary_charge / (crystal_temperature * SolidStateDetectors.kB), zs, rtol = 0.01))
end