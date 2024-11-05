# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors
using SolidStateDetectors: getVe, getVh, Vl, get_path_to_example_config_files
using StaticArrays
using LinearAlgebra
using Unitful

include("test_utils.jl")

T = Float32

sim = Simulation{T}(SSD_examples[:InvertedCoax])
timed_simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1, 0.05], verbose = false)

pos = CartesianPoint{T}(0.02,0,0.05); Edep = 1u"eV"
nbcc = NBodyChargeCloud(pos, Edep, 40, radius = T(0.0005), number_of_shells = 2)

@timed_testset "Diffusion and Self-Repulsion" begin
    evt = Event(nbcc)
    timed_simulate!(evt, sim, self_repulsion = true, diffusion = true)
    signalsum = T(0)
    for i in 1:length(evt.waveforms)
        signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
    end
    signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
    @info signalsum
    @test isapprox( signalsum, T(2), atol = 5e-3 )
end

@timed_testset "Charge Trapping" begin
    sim.detector = SolidStateDetector(sim.detector, BoggsChargeTrappingModel{T}())
    evt = Event(pos, Edep)
    timed_simulate!(evt, sim)
    signalsum = T(0)
    for i in 1:length(evt.waveforms)
        signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
    end
    signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
    @info signalsum
    @test signalsum < T(1.99)
end

@timed_testset "Test ADLChargeDriftModel" begin

    config_dict = Dict(
        "High-purity germanium" => joinpath(get_path_to_example_config_files(), "ADLChargeDriftModel/drift_velocity_config.yaml"),
        "Silicon" => joinpath(get_path_to_example_config_files(), "ADLChargeDriftModel/drift_velocity_Si_300K_config.yaml")
    )

    geom_sigdigits(::Type{Int})::Int = 12
    geom_sigdigits(::Type{Float32})::Int = 6
    geom_sigdigits(::Type{Float64})::Int = 12

    for T in (Float32, Float64) 
        @testset "Precision type: $(T)" begin
            
            for material in keys(config_dict)
                
                @testset "$material" begin
                    config = config_dict[material]
                    
                    @testset "x-axis aligned with <100> axis" begin

                        cdm = ADLChargeDriftModel(config, T=T, phi110 = T(π/4)); # <100> aligned with x-axis
                        E = collect(T, 10 .^(0:0.1:6));

                        r100 = @SVector T[1,0,0]
                        r111 = @SVector T[sqrt(1/3), sqrt(1/3), sqrt(1/3)]

                        @testset "Electrons <100>" begin
                            ve100_cdm = broadcast(E -> abs(getVe(E * r100, cdm) ⋅ r100), E)
                            ve100_true = map(E -> Vl(E,cdm.electrons.axis100), E)
                            @test isapprox(ve100_cdm, ve100_true, rtol = 10.0^(-geom_sigdigits(T)))
                        end

                        @testset "Electrons <111>" begin
                            ve111_cdm = broadcast(E -> abs(getVe(E * r111, cdm) ⋅ r111), E)
                            ve111_true = map(E -> Vl(E,cdm.electrons.axis111), E)
                            @test isapprox(ve111_cdm, ve111_true, rtol = 10.0^(-geom_sigdigits(T)))
                        end

                        @testset "Holes <100>" begin
                            vh100_cdm = broadcast(E -> abs(getVh(E * r100, cdm) ⋅ r100), E)
                            vh100_true = map(E -> Vl(E,cdm.holes.axis100), E)
                            @test isapprox(vh100_cdm, vh100_true, rtol = 10.0^(-geom_sigdigits(T)))
                        end

                        @testset "Holes <111>" begin
                            vh111_cdm = broadcast(E -> abs(getVh(E * r111, cdm) ⋅ r111), E)
                            vh111_true = map(E -> Vl(E,cdm.holes.axis111), E)
                            @test isapprox(vh111_cdm, vh111_true, rtol = 10.0^(-geom_sigdigits(T)))
                        end
                    end
                    
                    @testset "x-axis and <100> axis intersect at an angle of -30°" begin

                        cdm = ADLChargeDriftModel(config, T=T, phi110 = T(π/12)); # <110> and y-axis intersect at an angle of 15°
                        E = collect(T, 10 .^(0:0.1:6));

                        r100 = @SVector T[sqrt(3/4),-1/2,0]
                        r111 = @SVector T[(1 + sqrt(1/3))/2, (1 - sqrt(1/3))/2, sqrt(1/3)]

                        @testset "Electrons <100>" begin
                            ve100_cdm = broadcast(E -> abs(getVe(E * r100, cdm) ⋅ r100), E)
                            ve100_true = map(E -> Vl(E,cdm.electrons.axis100), E)
                            @test isapprox(ve100_cdm, ve100_true, rtol = 10.0^(-geom_sigdigits(T)))
                        end

                        @testset "Electrons <111>" begin
                            ve111_cdm = broadcast(E -> abs(getVe(E * r111, cdm) ⋅ r111), E)
                            ve111_true = map(E -> Vl(E,cdm.electrons.axis111), E)
                            @test isapprox(ve111_cdm, ve111_true, rtol = 10.0^(-geom_sigdigits(T)))
                        end

                        @testset "Holes <100>" begin
                            vh100_cdm = broadcast(E -> abs(getVh(E * r100, cdm) ⋅ r100), E)
                            vh100_true = map(E -> Vl(E,cdm.holes.axis100), E)
                            @test isapprox(vh100_cdm, vh100_true, rtol = 10.0^(-geom_sigdigits(T)))
                        end

                        @testset "Holes <111>" begin
                            vh111_cdm = broadcast(E -> abs(getVh(E * r111, cdm) ⋅ r111), E)
                            vh111_true = map(E -> Vl(E,cdm.holes.axis111), E)
                            @test isapprox(vh111_cdm, vh111_true, rtol = 10.0^(-geom_sigdigits(T)))
                        end
                    end
                end
            end
        
            config = joinpath(@__DIR__,"../examples/example_config_files/ADLChargeDriftModel/drift_velocity_config_axes.yaml")

            @testset "All axes in one plane" begin

                cdm = ADLChargeDriftModel(config, T=T); # <110> and y-axis intersect at an angle of 15°
                E = collect(T, 10 .^(0:0.1:6));

                r100 = @SVector T[1,0,0]
                r111 = @SVector T[sqrt(1/3),sqrt(2/3),0]

                @testset "Electrons <100>" begin
                    ve100_cdm = broadcast(E -> abs(getVe(E * r100, cdm) ⋅ r100), E)
                    ve100_true = map(E -> Vl(E,cdm.electrons.axis100), E)
                    @test isapprox(ve100_cdm, ve100_true, rtol = 10.0^(-geom_sigdigits(T)))
                end

                @testset "Electrons <111>" begin
                    ve111_cdm = broadcast(E -> abs(getVe(E * r111, cdm) ⋅ r111), E)
                    ve111_true = map(E -> Vl(E,cdm.electrons.axis111), E)
                    @test isapprox(ve111_cdm, ve111_true, rtol = 10.0^(-geom_sigdigits(T)))
                end

                @testset "Holes <100>" begin
                    vh100_cdm = broadcast(E -> abs(getVh(E * r100, cdm) ⋅ r100), E)
                    vh100_true = map(E -> Vl(E,cdm.holes.axis100), E)
                    @test isapprox(vh100_cdm, vh100_true, rtol = 10.0^(-geom_sigdigits(T)))
                end

                @testset "Holes <111>" begin
                    vh111_cdm = broadcast(E -> abs(getVh(E * r111, cdm) ⋅ r111), E)
                    vh111_true = map(E -> Vl(E,cdm.holes.axis111), E)
                    @test isapprox(vh111_cdm, vh111_true, rtol = 10.0^(-geom_sigdigits(T)))
                end
            end
        end
    end

    @testset "Test parsing of ADLChargeDriftModel config files with units" begin
        cdm0 = ADLChargeDriftModel() # default charge drift model
        @test cdm0.electrons.axis100.mu0  == 3.8609f0
        @test cdm0.electrons.axis100.beta == 0.805f0
        @test cdm0.electrons.axis100.E0   == 51100f0
        @test cdm0.electrons.axis100.mun  == -0.0171f0
        @test cdm0.electrons.axis111.mu0  == 3.8536f0
        @test cdm0.electrons.axis111.beta == 0.641f0
        @test cdm0.electrons.axis111.E0   == 53800f0
        @test cdm0.electrons.axis111.mun  == 0.051f0
        @test cdm0.holes.axis100.mu0  == 6.1824f0
        @test cdm0.holes.axis100.beta == 0.942f0
        @test cdm0.holes.axis100.E0   == 18500f0
        @test cdm0.holes.axis111.mu0  == 6.1215f0
        @test cdm0.holes.axis111.beta == 0.662f0
        @test cdm0.holes.axis111.E0   == 18200f0

        cdm_nounits = ADLChargeDriftModel(joinpath(get_path_to_example_config_files(), "ADLChargeDriftModel/drift_velocity_config_nounits.yaml"))
        @test cdm0.electrons == cdm_nounits.electrons
        @test cdm0.holes == cdm_nounits.holes
        @test cdm0.crystal_orientation ≈ cdm_nounits.crystal_orientation

    end

    @testset "Modify mobility parameters using keyword arguments" begin
        cdm = ADLChargeDriftModel{Float64}(e100μ0 = 40000u"cm^2/(V*s)", e111μn = 5u"cm^2/(V*s)", h100β = 2, h111E0 = 50u"V/cm")
        @test cdm.electrons.axis100.mu0 == 4.0     # internal unit is m^2/(V*s)
        @test cdm.electrons.axis111.mun == 0.0005  # internal unit is m^2/(V*s)
        @test cdm.holes.axis100.beta == 2.0        # no unit conversion
        @test cdm.holes.axis111.E0 == 5000.0       # internal unit is V/m
    end

    @testset "Construct ADLChargeDriftModel from dictionary (with units)" begin
        config = Dict(
            "phi110" => -45u"°",
            "material" => "HPGe",
            "drift" => Dict(
                "velocity" => Dict(
                    "parameters" => Dict(
                        "e100" => Dict(
                            "mu0" => 38609u"cm^2/(V*s)",
                            "beta" => 0.805,
                            "E0" => 511u"V/cm",
                            "mun" => -171u"cm^2/(V*s)"
                        ),
                        "e111" => Dict(
                            "mu0" => 38536u"cm^2/(V*s)",
                            "beta" => 0.641,
                            "E0" => 538u"V/cm",
                            "mun" => 510u"cm^2/(V*s)"
                        ),
                        "h100" => Dict(
                            "mu0" => 61824u"cm^2/(V*s)",
                            "beta" => 0.942,
                            "E0" => 185u"V/cm"
                        ),
                        "h111" => Dict(
                            "mu0" => 61215u"cm^2/(V*s)",
                            "beta" => 0.662,
                            "E0" => 182u"V/cm"
                        )
                    )
                )
            )
        )

        cdm0 = ADLChargeDriftModel() # default charge drift model
        cdmdict = ADLChargeDriftModel(config)

        @test cdm0.electrons == cdmdict.electrons
        @test cdm0.holes == cdmdict.holes
        @test cdm0.crystal_orientation == cdmdict.crystal_orientation
    end
end