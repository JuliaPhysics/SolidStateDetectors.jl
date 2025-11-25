using SolidStateDetectors
using Test

using Unitful

T = Float32

@timed_testset "Test impurity densities" begin
    @testset "Constant impurity density" begin
        d = Dict("impurity_density" => Dict(
                "name" => "constant",
                "value" => -5e9u"cm^-3"
            )
        )
        cd = SolidStateDetectors.ImpurityDensity(T, d["impurity_density"], SolidStateDetectors.default_unit_tuple())
        @test cd isa ConstantImpurityDensity{T}
        @test cd.ρ == -5f15
        @test cd == ConstantImpurityDensity(-5f9u"cm^-3")

        # passing an incompatible unit will throw a ConfigFileError
        @test_throws SolidStateDetectors.ConfigFileError ConstantImpurityDensity{T}(-5u"K")
    end 
    @testset "Linear impurity density" begin 
        d = Dict("impurity_density" => Dict(
                "name" => "linear",
                    "x" => Dict(
                        "init" => 1e-10,
                        "gradient" => 1.0e-11
                )
            )
        )
        cd = SolidStateDetectors.ImpurityDensity(T, d["impurity_density"], SolidStateDetectors.default_unit_tuple())
        @test cd isa LinearImpurityDensity{T}
        @test cd.offsets[1] == 1f-10
        @test cd.gradients[1] == 1f-11
    end
    @testset "Cylindrical impurity density" begin 
        d = Dict("impurity_density" => Dict(
                "name" => "cylindrical",
                    "r" => Dict(
                        "init" => 1e-10,
                        "gradient" => 1.0e-11
                )
            )
        )
        cd = SolidStateDetectors.ImpurityDensity(T, d["impurity_density"], SolidStateDetectors.default_unit_tuple())
        @test cd isa SolidStateDetectors.CylindricalImpurityDensity{T}
        @test cd.offsets[1] == 1f-10
        @test cd.gradients[1] == 1f-11
    end
end

@timed_testset "Test charge densities" begin
    @testset "Linear charge density" begin 
        d = Dict("charge_density" => Dict(
                "name" => "linear",
                    "x" => Dict(
                        "init" => 1e-10,
                        "gradient" => 1.0e-11
                )
            )
        )
        cd = SolidStateDetectors.ChargeDensity(T, d["charge_density"], SolidStateDetectors.default_unit_tuple())
        @test cd isa SolidStateDetectors.LinearChargeDensity{T}
        @test cd.offsets[1] == 1f-10
        @test cd.gradients[1] == 1f-11
    end
    @testset "Cylindrical charge density" begin 
        d = Dict("charge_density" => Dict(
                "name" => "cylindrical",
                    "r" => Dict(
                        "init" => 1e-10,
                        "gradient" => 1.0e-11
                )
            )
        )
        cd = SolidStateDetectors.ChargeDensity(T, d["charge_density"], SolidStateDetectors.default_unit_tuple())
        @test cd isa SolidStateDetectors.CylindricalChargeDensity{T}
        @test cd.offsets[1] == 1f-10
        @test cd.gradients[1] == 1f-11
    end
end

@timed_testset "Test boule impurity densities and corrections" begin
    sim = Simulation{T}("test_config_files/BEGe_01.yaml")

    @test SolidStateDetectors.get_impurity_density(sim.detector.semiconductor.impurity_density_model, CylindricalPoint{T}(0,0,0)) == T(0.8*T(SolidStateDetectors.to_internal_units(-5e9u"cm^-3")) + T(SolidStateDetectors.to_internal_units(-5e8u"cm^-3")))

    det_z0 = T(0.12)
    boule_ρ0 = T(-1e16)
    boule_gradient = T(-1e17)
    idm = LinBouleImpurityDensity{T}(boule_ρ0, boule_gradient, det_z0)
    sim.detector = SolidStateDetector(sim.detector, idm)

    det_ρ0 = SolidStateDetectors.get_impurity_density(sim.detector.semiconductor.impurity_density_model, CartesianPoint{T}(0,0,0))

    @test det_ρ0 == boule_ρ0 + boule_gradient * det_z0

    sim2 = Simulation{T}("test_config_files/BEGe_02.yaml")

    @test sim.detector.semiconductor.impurity_density_model == sim2.detector.semiconductor.impurity_density_model

    det_z0 = T(0.05)
    boule_ρ0 = T(-3e15)
    boule_gradient = T(-1e16)
    boule_n = T(-2e15)
    boule_l = T(0.05)
    boule_m = T(0.01)

    zimp = T.(collect(0:0.01:0.06))
    yimp = T.([ -3.013e15, -3.137e15, -3.3e15, -3.571e15, -4.136e15, -5.5e15, -9.037e15])

    idm = LinExpBouleImpurityDensity{T}(boule_ρ0, boule_gradient, boule_n, boule_l, boule_m, det_z0)
    sim.detector = SolidStateDetector(sim.detector, idm)
    
    det_ρ0 = SolidStateDetectors.get_impurity_density(sim.detector.semiconductor.impurity_density_model, CartesianPoint{T}(0,0,0))
    
    @test det_ρ0 == boule_ρ0 + boule_gradient * det_z0 + boule_n * exp((det_z0 - boule_l)/boule_m)

    timed_calculate_electric_potential!(sim, refinement_limits=0.01)
    U_est = timed_estimate_depletion_voltage(sim)

    idm_spline = SplineBouleImpurityDensity{T}(zimp, yimp, det_z0)
    sim.detector = SolidStateDetector(sim.detector, idm_spline)
    timed_calculate_electric_potential!(sim, refinement_limits=0.01)
    U_est_spline = timed_estimate_depletion_voltage(sim)

    @test isapprox(U_est, U_est_spline; atol=10u"V")
end
