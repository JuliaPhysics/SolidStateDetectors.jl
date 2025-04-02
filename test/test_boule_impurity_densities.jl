using SolidStateDetectors
using Test

T = Float32

@timed_testset "Test boule impurity densities" begin
    sim = Simulation{T}("BEGe_01.yaml")
    det_z0 = T(0.12)
    boule_ρ0 = T(-1e16)
    boule_gradient = T(-1e17)
    idm = LinBouleImpurityDensity{T}(boule_ρ0, boule_gradient, det_z0)
    sim.detector = SolidStateDetector(sim.detector, idm)

    det_ρ0 = SolidStateDetectors.get_impurity_density(sim.detector.semiconductor.impurity_density_model, CylindricalPoint{T}(0,0,0))

    @test det_ρ0 == boule_ρ0 + boule_gradient * det_z0

    boule_c = T(-2e15)
    boule_L = T(0.05)
    boule_tau = T(0.03)
    idm = LinExpBouleImpurityDensity{T}(boule_ρ0, boule_gradient, boule_c, boule_L, boule_tau, det_z0)
    sim.detector = SolidStateDetector(sim.detector, idm)

    det_ρ0 = SolidStateDetectors.get_impurity_density(sim.detector.semiconductor.impurity_density_model, CylindricalPoint{T}(0,0,0))
    
    @test det_ρ0 == boule_ρ0 + boule_gradient * det_z0 + boule_c * exp((det_z0 - boule_L)/boule_tau)
end
