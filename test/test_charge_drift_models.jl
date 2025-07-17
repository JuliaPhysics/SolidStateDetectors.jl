# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors
using SolidStateDetectors: getVe, getVh, Vl, get_path_to_example_config_files, AbstractChargeDriftModel, group_points_by_distance, distance_squared
using ArraysOfArrays
using InteractiveUtils
using StaticArrays
using LinearAlgebra
using Unitful

# include("test_utils.jl")

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

    config_dict = SolidStateDetectors.parse_config_file(SSD_examples[:InvertedCoax])
    @testset "Parse config file 1" begin
        config_dict["detectors"][1]["semiconductor"]["charge_trapping_model"] = Dict(
            "model" => "Boggs",
            "parameters" => Dict(
                "nσe" => "0.001cm^-1",
                "nσh" => "0.0005cm^-1",
                "temperature" => "78K"
            )
        )
        simA = @test_nowarn Simulation{T}(config_dict)
        @test simA.detector.semiconductor.charge_trapping_model isa BoggsChargeTrappingModel{T}
        @test simA.detector.semiconductor.charge_trapping_model.nσe == T(0.1)
        @test simA.detector.semiconductor.charge_trapping_model.nσh == T(0.05)
        @test simA.detector.semiconductor.charge_trapping_model.temperature == T(78)
    end
    @testset "Parse config file 2" begin
        config_dict["detectors"][1]["semiconductor"]["charge_trapping_model"] = Dict(
            "model" => "Boggs",
            "parameters" => Dict(
                "nσe-1" => "500cm",
                "nσh-1" => "500cm",
                "meffe" => 0.1,
                "meffh" => 0.2,
                "temperature" => "100K"
            )
        )
        simB = @test_nowarn Simulation{T}(config_dict)
        @test simB.detector.semiconductor.charge_trapping_model isa BoggsChargeTrappingModel{T}
        @test simB.detector.semiconductor.charge_trapping_model.nσe == T(0.2)
        @test simB.detector.semiconductor.charge_trapping_model.nσh == T(0.2)
        @test simB.detector.semiconductor.charge_trapping_model.meffe == T(0.1)
        @test simB.detector.semiconductor.charge_trapping_model.meffh == T(0.2)
        @test simB.detector.semiconductor.charge_trapping_model.temperature == T(100)

    end
end

@timed_testset "Test completeness of charge drift models" begin
    for c in InteractiveUtils.subtypes(AbstractChargeDriftModel)
        if isstructtype(c)
            @testset "$(c)" begin 
                cdm = c{T}
                fv = @SVector T[1,0,0]

                # default constructor
                @test cdm() isa c
                # electron-drift velocity
                @test hasmethod(getVe, Tuple{SVector{3,T}, cdm})
                @test getVe(fv, cdm()) isa SVector
                # hole-drift velocity
                @test hasmethod(getVh, Tuple{SVector{3,T}, cdm})
                @test getVh(fv, cdm()) isa SVector
            end
        end
    end
end


@timed_testset "Test IsotropicChargeDriftModel" begin

    @testset "Test constructors of IsotropicChargeDriftModel" begin
        cdm0 = IsotropicChargeDriftModel{T}() # default charge drift model
        @test cdm0.μ_e  == 0.1f0
        @test cdm0.μ_h == 0.1f0

        cdm1 = IsotropicChargeDriftModel{T}(1000u"cm^2/(V*s)", 1000u"cm^2/(V*s)")
        @test cdm1 == cdm0

        cdm2 = IsotropicChargeDriftModel{T}(0.1, 0.1)
        @test cdm2 == cdm0

        config_dict = Dict(
            "model" => "IsotropicChargeDriftModel",
            "mobilities" => Dict(
                "e" => "1000cm^2/(V*s)",
                "h" => "1000cm^2/(V*s)"
            )
        )
        cdm3 = IsotropicChargeDriftModel{T}(config_dict)
        @test cdm3 == cdm0
    end
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

    @testset "Test parsing of ADL2016ChargeDriftModel config files with units" begin
        cdm0 = ADL2016ChargeDriftModel(T=T) # default charge drift model
        @test cdm0.electrons.mu0      == 3.7165f0
        @test cdm0.electrons.beta     == 0.804f0
        @test cdm0.electrons.E0       == 50770f0
        @test cdm0.electrons.mun      == -0.0145f0
        @test cdm0.parameters.Γ0      == 0.496f0  # η0
        @test cdm0.parameters.Γ1      == 0.0296f0 # b
        @test cdm0.parameters.Γ2      == 120000f0 # Eref
        @test cdm0.holes.axis100.mu0  == 6.2934f0
        @test cdm0.holes.axis100.beta == 0.735f0
        @test cdm0.holes.axis100.E0   == 18190f0
        @test cdm0.holes.axis111.mu0  == 6.2383f0
        @test cdm0.holes.axis111.beta == 0.749f0
        @test cdm0.holes.axis111.E0   == 14390f0
    end

    @testset "Test equivalence of longitudinal drift parameter implementation" begin
    
        # The function to determine the hole drift for both models is equivalent,
        # so the hole drift parameters should be stored in the same way
        cdm2016 = ADL2016ChargeDriftModel(T=T)
        cdm = ADLChargeDriftModel(T=T,
            e100μ0 = 37165u"cm^2/(V*s)",
            e100β  = 0.804,
            e100E0 = 507.7u"V/cm",
            e100μn = -145u"cm^2/(V*s)",
            h100μ0 = 62934u"cm^2/(V*s)",
            h100β  = 0.735,
            h100E0 = 181.9u"V/cm",
            h111μ0 = 62383u"cm^2/(V*s)",
            h111β  = 0.749,
            h111E0 = 143.9u"V/cm"
        )
        
        @test cdm2016.electrons == cdm.electrons.axis100
        @test cdm2016.holes     == cdm.holes
    end
end

@timed_testset "Test grouping of charges" begin
    for N in 1:100
        @timed_testset "Test for length $(N)" begin
            
            # Generate random data
            pts = [CartesianPoint{T}(rand(3)...) for _ in Base.OneTo(N)]
            E = rand(T, N)
            d = rand(T)
            
            # Evaluate method
            ptsg, Eg = group_points_by_distance(pts, E, d)
            
            # Test correctness
            s0 = Set(pts)

            # result should be a VectorOfVectors
            @test ptsg isa VectorOfVectors{<:CartesianPoint{T}}

            # all points should appear in the result
            @test Set(flatview(ptsg))    == Set(pts)
            @test length(flatview(ptsg)) == length(pts)
            @test Set(flatview(Eg))      == Set(E)
            @test length(flatview(Eg))   == length(E)
        
            for (i,group) in enumerate(ptsg)
                @testset "Length $(N), subgroup $(i)" begin
                    sg = Set(group)
                    sc = setdiff(s0, group)
                    for pt in group
                        swo = setdiff(sg, Set([pt]))
                        @test isempty(swo) || any(distance_squared.(Ref(pt), swo) .<= d^2)
                        @test all(distance_squared.(Ref(pt), sc) .> d^2)
                    end
                end
            end
        end
    end
end