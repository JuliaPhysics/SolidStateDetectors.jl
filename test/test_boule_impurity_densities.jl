using SolidStateDetectors
using Test

T = Float32

@timed_testset "Test boule impurity densities and corrections" begin
    sim = Simulation{T}("BEGe_01.yaml")

    @test SolidStateDetectors.get_impurity_density(sim.detector.semiconductor.impurity_density_model, CylindricalPoint{T}(0,0,0)) == T(0.8*T(SolidStateDetectors.to_internal_units(-5e9u"cm^-3")) + T(SolidStateDetectors.to_internal_units(-5e8u"cm^-3")))

    det_z0 = T(0.12)
    boule_ρ0 = T(-1e16)
    boule_gradient = T(-1e17)
    idm = LinBouleImpurityDensity{T}(boule_ρ0, boule_gradient, det_z0)
    sim.detector = SolidStateDetector(sim.detector, idm)

    det_ρ0 = SolidStateDetectors.get_impurity_density(sim.detector.semiconductor.impurity_density_model, CylindricalPoint{T}(0,0,0))

    @test det_ρ0 == boule_ρ0 + boule_gradient * det_z0

    sim2 = Simulation{T}("BEGe_02.yaml")

    @test sim.detector.semiconductor.impurity_density_model == sim2.detector.semiconductor.impurity_density_model

    boule_n = T(-2e15)
    boule_l = T(0.05)
    boule_m = T(0.03)
    idm = LinExpBouleImpurityDensity{T}(boule_ρ0, boule_gradient, boule_n, boule_l, boule_m, det_z0)
    sim.detector = SolidStateDetector(sim.detector, idm)

    det_ρ0 = SolidStateDetectors.get_impurity_density(sim.detector.semiconductor.impurity_density_model, CylindricalPoint{T}(0,0,0))
    
    @test det_ρ0 == boule_ρ0 + boule_gradient * det_z0 + boule_n * exp((det_z0 - boule_l)/boule_m)
end
