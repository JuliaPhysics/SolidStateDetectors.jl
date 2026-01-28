using SolidStateDetectors
using Test
using StaticArrays
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
        pt = CartesianPoint{T}(1.0, 0.0, 0.0)
        @test SolidStateDetectors.get_impurity_density(cd, pt) == -5f15
        @test cd == ConstantImpurityDensity(-5f9u"cm^-3")
        cd_scaled = @test_nowarn (cd * 1.2)
        @test cd_scaled.ρ ≈ -6f15
        cd_offset = @test_nowarn cd + 1f9u"cm^-3"
        @test cd_offset.ρ ≈ -4f15

        # passing an incompatible unit will throw a ConfigFileError
        @test_throws SolidStateDetectors.ConfigFileError ConstantImpurityDensity{T}(-5u"K")
    end 
    @testset "Linear impurity density" begin 
        d = Dict("impurity_density" => Dict(
                "name" => "linear",
                "offset" => 1e10,
                "gradient" => Dict(
                    "x" => 1.0e11
                )
            )
        )
        cd = SolidStateDetectors.ImpurityDensity(T, d["impurity_density"], SolidStateDetectors.default_unit_tuple())
        @test cd isa LinearImpurityDensity{T}
        @test cd.offset == 1f10
        @test cd.gradients[1] == 1f11
        cd_scaled = @test_nowarn (cd * 1.2)
        @test cd_scaled.offset ≈ 1.2f10
        @test cd_scaled.gradients[1] ≈ 1.2f11
        cd_offset = @test_nowarn cd + 1f10u"m^-3"
        @test cd_offset.offset ≈ 2f10
        @test cd_offset.gradients == cd.gradients

        # should throw a warning, but still work and give the same results
        d_deprecated = Dict("impurity_density" => Dict(
                "name" => "linear",
                "x" => Dict(
                    "init" => 1e10,
                    "gradient" => 1.0e11
                )
            )
        )
        cd_deprecated = @test_logs (:warn,) SolidStateDetectors.ImpurityDensity(T, d_deprecated["impurity_density"], SolidStateDetectors.default_unit_tuple())
        @test cd_deprecated == cd

        pt = CartesianPoint{T}(1.0, 0.0, 0.0)
        ρ = SolidStateDetectors.get_impurity_density(cd, pt)
        # Expected: ρ = offset_x + x*grad_x = 1e10 + 1e11
        @test ρ ≈ T(1.1e11)

        # Impurity densities with not gradient should be updated to constant
        dc_deprecated = Dict("impurity_density" => Dict(
                "name" => "linear",
                "y" => Dict(
                    "init" => 1e10,
                    "gradient" => 0
                )
            )
        )

        cdc_deprecated = @test_logs (:warn,) SolidStateDetectors.ImpurityDensity(T, dc_deprecated["impurity_density"], SolidStateDetectors.default_unit_tuple())
        @test cdc_deprecated == ConstantImpurityDensity{T}(1e10)
    end
    @testset "Cylindrical impurity density" begin 
        d = Dict("impurity_density" => Dict(
                "name" => "cylindrical",
                "offset" => 1e10,
                "gradient" => Dict(
                    "r" => 1.0e11
                )
            )
        )
        cd = SolidStateDetectors.ImpurityDensity(T, d["impurity_density"], SolidStateDetectors.default_unit_tuple())
        @test cd isa SolidStateDetectors.CylindricalImpurityDensity{T}
        @test cd.offset == 1f10
        @test cd.gradients[1] == 1f11
        cd_scaled = @test_nowarn (cd * 1.2)
        @test cd_scaled.offset ≈ 1.2f10
        @test cd_scaled.gradients[1] ≈ 1.2f11
        cd_offset = @test_nowarn cd + 1f10u"m^-3"
        @test cd_offset.offset ≈ 2f10
        @test cd_offset.gradients == cd.gradients

        d_deprecated = Dict("impurity_density" => Dict(
                "name" => "cylindrical",
                "r" => Dict(
                    "init" => 1e10,
                    "gradient" => 1.0e11
                )
            )
        )
        cd_deprecated = @test_logs (:warn,) SolidStateDetectors.ImpurityDensity(T, d_deprecated["impurity_density"], SolidStateDetectors.default_unit_tuple())
        @test cd_deprecated == cd

        pt = CartesianPoint{T}(1.0, 0.0, 0.0)
        ρ = SolidStateDetectors.get_impurity_density(cd, pt)
        # Expected: ρ = offset_x + x*grad_x = 1e10 + 1e11
        @test ρ ≈ T(1.1e11)
    end
end

@timed_testset "Test charge densities" begin
    @testset "Linear charge density" begin 
        d = Dict("charge_density" => Dict(
                "name" => "linear",
                "offset" => 1e-10,
                "gradient" => Dict(
                    "x" => 1.0e-11
                )
            )
        )
        cd = SolidStateDetectors.ChargeDensity(T, d["charge_density"], SolidStateDetectors.default_unit_tuple())
        @test cd isa SolidStateDetectors.LinearChargeDensity{T}
        @test cd.offset == 1f-10
        @test cd.gradients[1] == 1f-11

        d_deprecated = Dict("charge_density" => Dict(
                "name" => "linear",
                "x" => Dict(
                    "init" => 1e-10,
                    "gradient" => 1.0e-11
                )
            )
        )
        cd_deprecated = @test_logs (:warn,) SolidStateDetectors.ChargeDensity(T, d_deprecated["charge_density"], SolidStateDetectors.default_unit_tuple())
        @test cd_deprecated == cd

        pt = CartesianPoint{T}(1.0, 0.0, 0.0)  
        ρ = SolidStateDetectors.get_charge_density(cd, pt)
        # Expected: ρ = offset_x + x*grad_x = 1e-10 + 1e-11
        @test ρ ≈ T(1.1e-10)
    end
    @testset "Cylindrical charge density" begin 
        d = Dict("charge_density" => Dict(
                "name" => "cylindrical",
                "offset" => 1e-10,
                "gradient" => Dict(
                    "r" => 1.0e-11
                )
            )
        )
        cd = SolidStateDetectors.ChargeDensity(T, d["charge_density"], SolidStateDetectors.default_unit_tuple())
        @test cd isa SolidStateDetectors.CylindricalChargeDensity{T}
        @test cd.offset == 1f-10
        @test cd.gradients[1] == 1f-11

        d_deprecated = Dict("charge_density" => Dict(
                "name" => "cylindrical",
                "r" => Dict(
                    "init" => 1e-10,
                    "gradient" => 1.0e-11
                )
            )
        )
        cd_deprecated = @test_logs (:warn,) SolidStateDetectors.ChargeDensity(T, d_deprecated["charge_density"], SolidStateDetectors.default_unit_tuple())
        @test cd_deprecated == cd

        pt = CartesianPoint{T}(1.0, 0.0, 0.0)  
        ρ = SolidStateDetectors.get_charge_density(cd, pt)
        # Expected: ρ = offset_x + x*grad_x = 1e-10 + 1e-11
        @test ρ ≈ T(1.1e-10)
    end
end

@timed_testset "Boule impurity densities and corrections" begin

    sim = Simulation{T}("test_config_files/BEGe_01.yaml")
    @test SolidStateDetectors.get_impurity_density(sim.detector.semiconductor.impurity_density_model, CylindricalPoint{T}(0,0,0)) == T(0.8*T(SolidStateDetectors.to_internal_units(-5e9u"cm^-3")) + T(SolidStateDetectors.to_internal_units(-5e8u"cm^-3")))

    det_z0 = T(0.05)
    boule_ρ0 = T(-3e15)
    boule_gradient = T(-1e16)
    boule_c = T(-1e15)
    boule_n = T(-2e15)
    boule_l = T(0.05)
    boule_m = T(0.01)
    
    # LinBouleImpurityDensity
    d = Dict("impurity_density" => Dict(
        "name" => "linear_boule",
        "a" => boule_ρ0,
        "b" => boule_gradient,
        "det_z0" => det_z0
    ))

    idm = LinBouleImpurityDensity{T}(boule_ρ0, boule_gradient, det_z0)
    @test (1 * idm) + 0 == idm
    @test LinBouleImpurityDensity{T}([boule_ρ0, boule_gradient], det_z0) == idm
    @test SolidStateDetectors.ImpurityDensity(T, d["impurity_density"], SolidStateDetectors.default_unit_tuple()) == idm
    @test SolidStateDetectors.get_impurity_density(idm, CartesianPoint{T}(0,0,0)) == boule_ρ0 + boule_gradient * det_z0

    # LinExpBouleImpurityDensity
    d = Dict("impurity_density" => Dict(
        "name" => "linear_exponential_boule",
        "a" => boule_ρ0,
        "b" => boule_gradient,
        "n" => boule_n,
        "l" => boule_l,
        "m" => boule_m,
        "det_z0" => det_z0
    ))

    idm = LinExpBouleImpurityDensity{T}(boule_ρ0, boule_gradient, boule_n, boule_l, boule_m, det_z0)
    @test (1 * idm) + 0 == idm
    @test LinExpBouleImpurityDensity{T}([boule_ρ0, boule_gradient, boule_n, boule_l, boule_m], det_z0) == idm
    @test SolidStateDetectors.ImpurityDensity(T, d["impurity_density"], SolidStateDetectors.default_unit_tuple()) == idm
    @test SolidStateDetectors.get_impurity_density(idm, CartesianPoint{T}(0,0,0)) == boule_ρ0 + boule_gradient * det_z0 + boule_n * exp((det_z0 - boule_l)/boule_m)

    sim.detector = SolidStateDetector(sim.detector, idm)
    timed_calculate_electric_potential!(sim, refinement_limits=0.01)
    U_est = timed_estimate_depletion_voltage(sim)

    # ParBouleImpurityDensity
    d = Dict("impurity_density" => Dict(
        "name" => "parabolic_boule",
        "a" => boule_ρ0,
        "b" => boule_gradient,
        "c" => boule_c,
        "det_z0" => det_z0
    ))

    idm = ParBouleImpurityDensity{T}(boule_ρ0, boule_gradient, boule_c, det_z0)
    @test (1 * idm) + 0 == idm
    @test ParBouleImpurityDensity{T}([boule_ρ0, boule_gradient, boule_c], det_z0) == idm
    @test SolidStateDetectors.ImpurityDensity(T, d["impurity_density"], SolidStateDetectors.default_unit_tuple()) == idm
    @test SolidStateDetectors.get_impurity_density(idm, CartesianPoint{T}(0,0,0)) == boule_ρ0 + boule_gradient * det_z0 + boule_c * det_z0^2

    # ParExpBouleImpurityDensity
    d = Dict("impurity_density" => Dict(
        "name" => "parabolic_exponential_boule",
        "a" => boule_ρ0,
        "b" => boule_gradient,
        "c" => boule_c,
        "n" => boule_n,
        "l" => boule_l,
        "m" => boule_m,
        "det_z0" => det_z0
    ))

    idm = ParExpBouleImpurityDensity{T}(boule_ρ0, boule_gradient, boule_c, boule_n, boule_l, boule_m, det_z0)
    @test (1 * idm) + 0 == idm
    @test ParExpBouleImpurityDensity{T}([boule_ρ0, boule_gradient, boule_c, boule_n, boule_l, boule_m], det_z0) == idm
    @test SolidStateDetectors.ImpurityDensity(T, d["impurity_density"], SolidStateDetectors.default_unit_tuple()) == idm
    @test SolidStateDetectors.get_impurity_density(idm, CartesianPoint{T}(0,0,0)) == boule_ρ0 + boule_gradient * det_z0 + boule_c * det_z0^2 + boule_n * exp((det_z0 - boule_l)/boule_m)

    # SplineBouleImpurityDensity
    zimp = T.(collect(0:0.01:0.06))
    yimp = T.([ -3.013e15, -3.137e15, -3.3e15, -3.571e15, -4.136e15, -5.5e15, -9.037e15])
    idm_spline = SplineBouleImpurityDensity{T}(zimp, yimp, det_z0)
    @test ((1 * idm_spline) + 0).ρ == idm_spline.ρ
    sim.detector = SolidStateDetector(sim.detector, idm_spline)
    timed_calculate_electric_potential!(sim, refinement_limits=0.01)
    U_est_spline = timed_estimate_depletion_voltage(sim)
    @test isapprox(U_est, U_est_spline; atol=10u"V")
end

@testset "PtypePNJunctionImpurityDensity" begin
    d = Dict("impurity_density" => Dict(
        "name" => "PtypePNjunction",
        "lithium_annealing_temperature" => "623K",
        "lithium_annealing_time" => "18minute",
        "doped_contact_id" => 2,
        "bulk_impurity_density" => Dict(
            "name" => "constant",
            "value" => "-1e10cm^-3"
        )
    ))

    idm = SolidStateDetectors.ImpurityDensity(T, d["impurity_density"], SolidStateDetectors.default_unit_tuple())
    @test (1 * idm) + 0 == idm
    @test idm isa PtypePNJunctionImpurityDensity{T}
    @test idm.bulk_imp_model isa ConstantImpurityDensity{T}
    
    sidm = idm.surface_imp_model
    @test (1 * sidm) + 0 == sidm
    @test sidm isa ThermalDiffusionLithiumDensity{T}
    
    d["impurity_density"]["name"] = "li_diffusion"
    delete!(d["impurity_density"], "bulk_impurity_density")
    sidm2 = SolidStateDetectors.ImpurityDensity(T, d["impurity_density"], SolidStateDetectors.default_unit_tuple())
    @test sidm2 isa ThermalDiffusionLithiumDensity{T}
    @test sidm.lithium_density_on_contact == sidm2.lithium_density_on_contact
    @test sidm.lithium_diffusivity == sidm2.lithium_diffusivity
end

@testset "ThermalDiffusionLithiumParameters" begin
    # test all the error messages constructing a lithium density parameter dictionary
    d = Dict{String, Any}("model" => "ThermalDiffusionLithiumDensity")

    # missing annealing temperature ranges
    @test_throws SolidStateDetectors.ConfigFileError SolidStateDetectors.ThermalDiffusionLithiumParameters(d)

    # missing experimental parameters
    d["annealing_temperature_ranges"] = [
        Dict("T_min" => "100K", "T_max" => "200K", "D0" => "25e-4cm^2/s", "H" => "11800cal/mol"),
        Dict("T_min" => "200K", "T_max" => "400K", "D0" => "20e-4cm^2/s", "H" => "12300cal/mol")
    ]
    @test_throws SolidStateDetectors.ConfigFileError SolidStateDetectors.ThermalDiffusionLithiumParameters(d)

    # missing parameters a and b
    d["experimental_parameters"] = Dict{String,Any}()
    @test_throws SolidStateDetectors.ConfigFileError SolidStateDetectors.ThermalDiffusionLithiumParameters(d)
    d["experimental_parameters"]["a"] = 21.27
    @test_throws SolidStateDetectors.ConfigFileError SolidStateDetectors.ThermalDiffusionLithiumParameters(d)
    d["experimental_parameters"]["b"] = 2610
    lithium_parameters = @test_nowarn SolidStateDetectors.ThermalDiffusionLithiumParameters(d; T)
    @test lithium_parameters isa SolidStateDetectors.ThermalDiffusionLithiumDensityParameters{T}

    # calculate lithium diffusivity inside and outside of the temperature range
    @test SolidStateDetectors.calculate_lithium_diffusivity(T(200), lithium_parameters.diffusion) isa T
    @test_throws ArgumentError SolidStateDetectors.calculate_lithium_diffusivity(T(500), lithium_parameters.diffusion)

    # empty annealing temperature ranges
    d["annealing_temperature_ranges"] = []
    @test_throws SolidStateDetectors.ConfigFileError SolidStateDetectors.ThermalDiffusionLithiumParameters(d)

    # missing H
    d["annealing_temperature_ranges"] = [Dict("T_min" => "100K", "T_max" => "200K", "D0" => "25e-4cm^2/s")]
    @test_throws SolidStateDetectors.ConfigFileError SolidStateDetectors.ThermalDiffusionLithiumParameters(d)

    # extra key
    d["annealing_temperature_ranges"] = [Dict("T_min" => "100K", "T_max" => "200K", "D0" => "25e-4cm^2/s", "H" => "11800cal/mol", "extra" => "key")]
    @test_logs (:warn,) SolidStateDetectors.ThermalDiffusionLithiumParameters(d)

    # mixed up temperatures
    d["annealing_temperature_ranges"] = [Dict("T_min" => "200K", "T_max" => "100K", "D0" => "25e-4cm^2/s", "H" => "11800cal/mol")]
    @test_throws SolidStateDetectors.ConfigFileError SolidStateDetectors.ThermalDiffusionLithiumParameters(d)

    # missing temperatures
    d["annealing_temperature_ranges"] = [
        Dict("T_min" => "100K", "T_max" => "200K", "D0" => "25e-4cm^2/s", "H" => "11800cal/mol"),
        Dict("T_min" => "300K", "T_max" => "400K", "D0" => "20e-4cm^2/s", "H" => "12300cal/mol")
    ]
    @test_throws SolidStateDetectors.ConfigFileError SolidStateDetectors.ThermalDiffusionLithiumParameters(d)
end

@testset "Charge density types" begin
    # Test get_charge_density accepts AbstractCoordinatePoint
    lcd = SolidStateDetectors.LinearChargeDensity{Float32}(1f0, (0f0, 1f0, 0f0))
    ccd = SolidStateDetectors.CylindricalChargeDensity{Float32}(2f0, (1f0, 0f0, -1f0))

    cart_pt = CartesianPoint{Float32}(0.5f0, 0.0f0, 2f0)
    cyl_pt = CylindricalPoint{Float32}(0.5f0, 0.0f0, 1f0)

    @test @inferred(SolidStateDetectors.get_charge_density(lcd, cart_pt)) isa Float32
    @test @inferred(SolidStateDetectors.get_charge_density(ccd, cyl_pt)) isa Float32
end
