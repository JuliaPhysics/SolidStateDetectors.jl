# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors
using SolidStateDetectors: getVe, getVh, Vl, get_path_to_example_config_files, AbstractChargeDriftModel, ConstantImpurityDensity, group_points_by_distance, distance_squared, scale_to_temperature
using SolidStateDetectors.ConstructiveSolidGeometry: AbstractCoordinatePoint, geom_sigdigits
using SolidStateDetectors: AbstractVirtualVolume
using ArraysOfArrays
using InteractiveUtils
using StaticArrays
using LinearAlgebra
using Unitful

# include("test_utils.jl")

T = Float32

pos = CartesianPoint{T}(0.02,0,0.05); Edep = 1u"eV"
nbcc = NBodyChargeCloud(pos, Edep, 40, radius = T(0.0005), number_of_shells = 2)

@timed_testset "Create single events" begin
    evt0 = Event([pos], [Edep], 40, radius = [T(0.0005)], number_of_shells = 2)
    # Test out all alternative methods to create an Event
    for evt in (
        Event([pos], [Edep], 40, radius = [T(0.0005)], number_of_shells = 2, max_interaction_distance = 1u"m"),
        Event([[pos]], [[Edep]], 40, radius = [[T(0.0005)]], number_of_shells = 2),
        Event(nbcc),
        Event([nbcc]),
        Event([nbcc.locations], [nbcc.energies]),
        Event(nbcc.locations, nbcc.energies, max_interaction_distance = 1u"m")
    )
        @test evt.locations == evt0.locations
        @test evt.energies == evt0.energies
    end
end

sim = Simulation{T}(SSD_examples[:InvertedCoax])
timed_simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1, 0.05], verbose = false)

@timed_testset "Diffusion and Self-Repulsion" begin
    evt = Event(nbcc)
    timed_simulate!(evt, sim, self_repulsion = true, diffusion = true)
    signalsum = T(0)
    for i in 1:length(evt.waveforms)
        signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
    end
    signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
    @test isapprox( signalsum, T(2), atol = 5e-3 )
end

@timed_testset "Move charges inside semiconductor" begin

    # define a charge cloud very outside of the detector => should throw an AssertionError
    nbcc_top = NBodyChargeCloud(CartesianPoint{T}(0.01,0,0.081), Edep, 10, radius = T(0.002), number_of_shells = 1)
    evt = Event(nbcc_top)
    @test_throws AssertionError timed_simulate!(evt, sim, self_repulsion = true, diffusion = true)

    # define a charge cloud very close to the top surface of the example inverted coax detector (top = 80mm) => should throw a warning
    nbcc_top = NBodyChargeCloud(CartesianPoint{T}(0.01,0,0.079), Edep, 10, radius = T(0.002), number_of_shells = 1)
    evt = Event(nbcc_top)
    @test_logs (:warn, r".*Moving them inside.*") timed_simulate!(evt, sim, self_repulsion = true, diffusion = true)

    signalsum = T(0)
    for i in 1:length(evt.waveforms)
        signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
    end
    signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
    @test isapprox( signalsum, T(2), atol = 5e-3 )
end

@timed_testset "Initial radius for charge carriers" begin
    for ptype in InteractiveUtils.subtypes(SolidStateDetectors.ParticleType)
        @test SolidStateDetectors.radius_guess(T(1e6), ptype) isa T
    end
end

@timed_testset "Charge Trapping: BoggsChargeTrappingModel" begin
    sim.detector = SolidStateDetector(sim.detector, BoggsChargeTrappingModel{T}())
    evt = Event(pos, Edep)
    timed_simulate!(evt, sim)
    signalsum = T(0)
    for i in 1:length(evt.waveforms)
        signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
    end
    signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
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
@timed_testset "Charge Trapping: ConstantLifetimeChargeTrappingModel" begin
    # test 1: parse the lifetimes and the inactive layer geometry
    simA = @test_nowarn Simulation{T}(SSD_examples[:TrueCoaxial])
    simA.detector = SolidStateDetector(simA.detector, ConstantImpurityDensity{T}(-1e-10))
    @testset "Parse the TrueCoaxial config file" begin
        @test simA.detector.semiconductor.charge_trapping_model isa CombinedChargeTrappingModel{T}
        @test simA.detector.semiconductor.charge_trapping_model.bulk_charge_trapping_model.τh == T(1e-3)
        @test simA.detector.semiconductor.charge_trapping_model.bulk_charge_trapping_model.τe == T(1e-3)
        @test simA.detector.semiconductor.charge_trapping_model.inactive_charge_trapping_model.τh == T(1e-6)
        @test simA.detector.semiconductor.charge_trapping_model.inactive_charge_trapping_model.τe  == T(1e-6)
        @test simA.detector.semiconductor.charge_trapping_model.inactive_layer_geometry.origin == CartesianPoint{T}(0.0, 0.0, 0.005)
        r0, r1 = T.((0.008957282, 0.01))
        @test simA.detector.semiconductor.charge_trapping_model.inactive_layer_geometry.r == tuple((r0, r1), (r0, r1))
        @test simA.detector.semiconductor.charge_trapping_model.inactive_layer_geometry.hZ == T(0.005)
    end

    timed_simulate!(simA, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1, 0.05], verbose = false)
    simA_inactive_layer_geometry=deepcopy(simA.detector.semiconductor.charge_trapping_model.inactive_layer_geometry)

    # test2: only model the bulk volume, test the bulk signals while varying the lifetime: τ
    pos = CartesianPoint{T}(0.005,0,0.005); Edep = 1u"eV"
    signalsum_list = T[]
    τ_list = (10 .^ (-3:0.3:0))*u"ms"
    for τ in τ_list
        parameters = Dict("parameters" => Dict("τh" => τ, "τe" => τ))
        trapping_model=ConstantLifetimeChargeTrappingModel{T}(parameters)
        simA.detector = SolidStateDetector(simA.detector, trapping_model)
        evt = Event(pos, Edep)
        timed_simulate!(evt, simA)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(simA.detector.semiconductor.material)))
        push!(signalsum_list, signalsum)
    end
    @test all(signalsum_list .< T(2.0))
    @test all(diff(signalsum_list) .> 0)

    # test3: model both the bulk/inactive volumes
    ## test 3.1: the bulk/inactive signals while varying the lifetimes: τ, τ_inactive
    pos_bulk = CartesianPoint{T}(0.0085,0,0.005); Edep = 1u"eV"
    @test !in(pos_bulk, simA_inactive_layer_geometry)
    signalsum_list_bulk = T[]
    τ_list = (10 .^ (-2:0.2:0))*u"ms"
    pos_inactive = CartesianPoint{T}(0.0095,0,0.005); Edep = 1u"eV"
    @test in(pos_inactive, simA_inactive_layer_geometry)
    signalsum_list_inactive = T[]
    τ_inactive_list = τ_list/100
    for (τ,τ_inactive) in zip(τ_list, τ_inactive_list)
        parameters = Dict("model" => "ConstantLifetime", "model_inactive" => "ConstantLifetime", "parameters" => Dict("τh" => τ, "τe" => τ), "parameters_inactive" => Dict("τh" => τ_inactive, "τe" => τ_inactive),
            "inactive_layer_geometry" => simA_inactive_layer_geometry)
        trapping_model=CombinedChargeTrappingModel{T}(parameters)
        simA.detector = SolidStateDetector(simA.detector, trapping_model)
        evt_bulk = Event(pos_bulk , Edep)
        timed_simulate!(evt_bulk, simA)
        signalsum_bulk  = T(0)
        for i in 1:length(evt_bulk.waveforms)
            signalsum_bulk  += abs(ustrip(evt_bulk.waveforms[i].signal[end]))
        end
        signalsum_bulk  *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(simA.detector.semiconductor.material)))
        push!(signalsum_list_bulk, signalsum_bulk)

        evt_inactive = Event(pos_inactive , Edep)
        timed_simulate!(evt_inactive, simA)
        signalsum_inactive  = T(0)
        for i in 1:length(evt_inactive.waveforms)
            signalsum_inactive  += abs(ustrip(evt_inactive.waveforms[i].signal[end]))
        end
        signalsum_inactive  *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(simA.detector.semiconductor.material)))
        push!(signalsum_list_inactive, signalsum_inactive)
    end
    @test all(signalsum_list_bulk .< T(2.0))
    @test all(diff(signalsum_list_bulk) .> 0)
    @test all(signalsum_list_inactive .< T(2.0))
    @test all(diff(signalsum_list_inactive) .> 0)
    @test all(signalsum_list_bulk .> signalsum_list_inactive)

    ## test 3.2: the inactive layer signals while varying the depth
    τ, τ_inactive = 1u"ms", 100u"ns"
    parameters = Dict("model" => "ConstantLifetime", "model_inactive" => "ConstantLifetime", "parameters" => Dict("τh" => τ, "τe" => τ), "parameters_inactive" => Dict("τh" => τ_inactive, "τe" => τ_inactive),
        "inactive_layer_geometry" => simA_inactive_layer_geometry)
    trapping_model=CombinedChargeTrappingModel{T}(parameters)
    simA.detector = SolidStateDetector(simA.detector, trapping_model)
    signalsum_list_inactive = T[]
    for depth in (0.1:0.1:0.9)/1000
        pos_inactive = CartesianPoint{T}(0.01-depth,0,0.005); Edep = 1u"eV"
        @test in(pos_inactive, simA_inactive_layer_geometry)
        evt_inactive = Event(pos_inactive , Edep)
        timed_simulate!(evt_inactive, simA, Δt=1u"ns")
        signalsum_inactive  = T(0)
        for i in 1:length(evt_inactive.waveforms)
            signalsum_inactive  += abs(ustrip(evt_inactive.waveforms[i].signal[end]))
        end
        signalsum_inactive  *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(simA.detector.semiconductor.material)))
        push!(signalsum_list_inactive, signalsum_inactive)
    end
    @test all(signalsum_list_inactive .< T(2.0))
    @test all(diff(signalsum_list_inactive) .> 0)
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

@timed_testset "Test InactiveLayerChargeDriftModel" begin
    @testset "Test constructors of InactiveLayerChargeDriftModel" begin
        cdm0 = InactiveLayerChargeDriftModel{T}()
        @test cdm0 isa AbstractChargeDriftModel{T}
        @test cdm0.temperature == T(90)
        @test cdm0.neutral_imp_model == ConstantImpurityDensity{T}(1e21)
        @test cdm0.bulk_imp_model == ConstantImpurityDensity{T}(-1e16)
        @test cdm0.surface_imp_model == ConstantImpurityDensity{T}(0)
        # the charge drift model should define methods for coordinate points and charge carriers:
        @test hasmethod(SolidStateDetectors.calculate_mobility, Tuple{typeof(cdm0), CartesianPoint{T}, Type{SolidStateDetectors.Electron}})
        @test hasmethod(SolidStateDetectors.calculate_mobility, Tuple{typeof(cdm0), CartesianPoint{T}, Type{SolidStateDetectors.Hole}})
        @test SolidStateDetectors.calculate_mobility(cdm0, CartesianPoint{T}(0,0,0), SolidStateDetectors.Electron) isa T
        @test SolidStateDetectors.calculate_mobility(cdm0, CartesianPoint{T}(0,0,0), SolidStateDetectors.Hole) isa T
    end

    simA = @test_nowarn Simulation{T}(SSD_examples[:TrueCoaxial])
    @testset "Parse the TrueCoaxial config file" begin
        @test simA.detector.semiconductor.charge_drift_model isa AbstractChargeDriftModel{T}
        @test simA.detector.semiconductor.charge_drift_model.temperature == T(90)
        @test simA.detector.semiconductor.charge_drift_model.neutral_imp_model == ConstantImpurityDensity{T}(5.6769e21)
        @test simA.detector.semiconductor.charge_drift_model.bulk_imp_model == simA.detector.semiconductor.impurity_density_model.bulk_imp_model
        @test simA.detector.semiconductor.charge_drift_model.surface_imp_model == simA.detector.semiconductor.impurity_density_model.surface_imp_model
        # the charge drift model should define methods for coordinate points and charge carriers:
        @test hasmethod(SolidStateDetectors.calculate_mobility, Tuple{typeof(simA.detector.semiconductor.charge_drift_model), CartesianPoint{T}, Type{SolidStateDetectors.Electron}})
        @test hasmethod(SolidStateDetectors.calculate_mobility, Tuple{typeof(simA.detector.semiconductor.charge_drift_model), CartesianPoint{T}, Type{SolidStateDetectors.Hole}})
        @test SolidStateDetectors.calculate_mobility(simA.detector.semiconductor.charge_drift_model, CartesianPoint{T}(0,0,0), SolidStateDetectors.Electron) isa T
        @test SolidStateDetectors.calculate_mobility(simA.detector.semiconductor.charge_drift_model, CartesianPoint{T}(0,0,0), SolidStateDetectors.Hole) isa T
    end

    @testset "Test if the detector is depleted (inactive layer is not taken into account for depletion) and in_inactive_layer" begin
        mm = 1 / 1000
        pn_r = 8.957282 * mm
        g = Grid(simA)
        ax1, ax2, ax3 = g.axes
        bulk_tick_dis, dl_tick_dis  = 0.05 * mm, 0.01 * mm
        user_additional_ticks_ax1 = sort(vcat(ax1.interval.left:bulk_tick_dis:pn_r, pn_r:dl_tick_dis:ax1.interval.right))
        user_ax1 = typeof(ax1)(ax1.interval, user_additional_ticks_ax1)
        user_g = typeof(g)((user_ax1, ax2, ax3))
        calculate_electric_potential!(simA, grid=user_g, depletion_handling=true)
        @test is_depleted(simA.point_types)
        pt_active = CartesianPoint{T}(0.005, 0.0, 0.005)
        pt_inactive = CartesianPoint{T}(0.0095, 0.0, 0.005)
        @test SolidStateDetectors.in_inactive_layer(pt_active, nothing, simA.point_types) == false
        @test SolidStateDetectors.in_inactive_layer(pt_inactive, nothing, simA.point_types) == true
    end
end

@timed_testset "Test IsotropicChargeDriftModel" begin

    @testset "Test constructors of IsotropicChargeDriftModel" begin
        cdm0 = IsotropicChargeDriftModel{T}() # default charge drift model
        @test cdm0.μ_e  == 0.1f0
        @test cdm0.μ_h == 0.1f0

        cdm1 = IsotropicChargeDriftModel{T}(μ_e = 1000u"cm^2/(V*s)", μ_h = 1000u"cm^2/(V*s)")
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
        
            config = "test_config_files/drift_velocity_config_axes.yaml"

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
        cdm0 = ADLChargeDriftModel() # default charge drift model with units
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
        
        # default charge drift model config has units, input units should be ignored
        input_units = SolidStateDetectors.construct_units(Dict("units" => Dict("length" => "mm", "potential" => "V", "angle" => "deg", "temperature" => "K")))

        cdm0_inputunits = ADLChargeDriftModel(SolidStateDetectors.default_ADL_config_file, input_units)
        @test cdm0.electrons == cdm0_inputunits.electrons
        @test cdm0.holes == cdm0_inputunits.holes

        # charge drift model config has no units, internal units are assumed
        cdm_nounits = ADLChargeDriftModel("test_config_files/drift_velocity_config_nounits.yaml")
        @test cdm0.electrons == cdm_nounits.electrons
        @test cdm0.holes == cdm_nounits.holes
        @test cdm0.crystal_orientation ≈ cdm_nounits.crystal_orientation

        # charge drift model config has no units, input units will be applied
        cdm_nounits_inputunits = ADLChargeDriftModel("test_config_files/drift_velocity_config_nounits.yaml", input_units)
        @test cdm_nounits.electrons != cdm_nounits_inputunits.electrons
        @test cdm_nounits.holes != cdm_nounits_inputunits.holes
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

    @testset "Test temperature scaling of drift parameters" begin
    
        # If a temperature is given to the ADL2016 model, the ChargeDriftModel is scaled. But at 77K (the reference temperature),
        # the parameters should be identical to the unscaled ADL2016 model
        
        cdm_77 = ADL2016ChargeDriftModel(temperature = 77u"K")
        cdm_100 = ADL2016ChargeDriftModel(temperature = 100u"K")
        cdm = ADL2016ChargeDriftModel()
        
        @test cdm_77.electrons == cdm.electrons
        @test cdm_77.holes     == cdm.holes
        @test cdm_100.electrons.mu0 < cdm.electrons.mu0
        @test cdm_100.holes.axis100.mu0 < cdm.holes.axis100.mu0
        @test cdm_100.holes.axis111.mu0 < cdm.holes.axis111.mu0

        cdm_2scalings = scale_to_temperature(scale_to_temperature(cdm, 87u"K"), 100u"K")
        @test cdm_2scalings.temperaturemodel.reference_temperature == cdm_100.temperaturemodel.reference_temperature
        @test isapprox(cdm_2scalings.electrons.E0, cdm_100.electrons.E0, atol = 5e-2)
        @test isapprox(cdm_2scalings.electrons.mu0, cdm_100.electrons.mu0, atol = 5e-2)
    end
end

@timed_testset "Test grouping of charges" begin
    for N in 1:100
        @testset "Test for length $(N)" begin
            
            # Generate random data
            pts = [CartesianPoint{T}(rand(3)...) for _ in Base.OneTo(N)]
            E = rand(T, N)
            d = rand(T)
            
            # Evaluate method
            ptsg0    = group_points_by_distance(pts, d)
            ptsg, Eg = group_points_by_distance(pts, E, d)
            @test ptsg0 == ptsg
            
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

struct OutsideTestVolume{T} <: SolidStateDetectors.AbstractVirtualVolume{T} end
Base.in(::CartesianPoint{T}, ::OutsideTestVolume{T}) where {T} = false

@testset "Virtual volumes" begin
    struct DummyGeom end
    Base.in(::AbstractCoordinatePoint{T}, ::DummyGeom) where {T} = true
    struct DummyVV{T} <: AbstractVirtualVolume{T}
        geometry::DummyGeom
    end
    pt = CartesianPoint{Float64}(0, 0, 0)
    @test in(pt, DummyVV{Float64}(DummyGeom()))
end


@testset "_calculate_signal" begin
    # Test that the function accepts AbstractVector{CartesianPoint{T}} and AbstractVector{T}
    model = NoChargeTrappingModel{T}()
    path = [CartesianPoint{T}(0.0, 0.0, 0.0), CartesianPoint{T}(0.1, 0.0, 0.0)]
    times = T[0.0, 1.0]
    wpot = SolidStateDetectors.interpolated_scalarfield(sim.weighting_potentials[1])
    
    signal = SolidStateDetectors._calculate_signal(model, path, times, one(T), wpot, sim.point_types)
    @test signal isa Vector{T}
    @test length(signal) == length(times)
end

@testset "Modulate Drift Vector" begin
    T = Float64
    sv = CartesianVector{T}(1,2,3)
    pt = CartesianPoint{T}(0,0,0)

    example_primitive_dir = joinpath(@__DIR__, "../examples/example_primitive_files")
    geom = SolidStateDetectors.ConstructiveSolidGeometry.Geometry(T, joinpath(example_primitive_dir, "Box.yaml"))

    arb_vol = SolidStateDetectors.ArbitraryDriftModificationVolume{T}("arb", 1, geom)
    @test_throws ErrorException SolidStateDetectors.modulate_driftvector(sv, pt, [arb_vol])

    step_vectors = [CartesianVector{T}(1,0,0), CartesianVector{T}(0,1,0)]
    current_pos  = [pt, pt]
    @test_throws ErrorException SolidStateDetectors._modulate_driftvectors!(step_vectors, current_pos, [arb_vol])

    empty_vdv = SolidStateDetectors.AbstractVirtualVolume{T}[]
    @test SolidStateDetectors.modulate_driftvector(sv, pt, empty_vdv) == sv

    result = SolidStateDetectors.modulate_driftvector(sv, pt, [OutsideTestVolume{T}()])
    @test result == sv

    step_vectors = [CartesianVector{T}(1,0,0), CartesianVector{T}(0,1,0)]
    current_pos  = [pt, pt]
    SolidStateDetectors._modulate_driftvectors!(step_vectors, current_pos, [OutsideTestVolume{T}()])
    @test step_vectors == [CartesianVector{T}(1,0,0), CartesianVector{T}(0,1,0)]

    dead_vol = SolidStateDetectors.DeadVolume{T, typeof(geom)}("dead", geom)
    result = SolidStateDetectors.modulate_driftvector(sv, pt, dead_vol)
    @test result == CartesianVector{T}(0,0,0)

    result = SolidStateDetectors.modulate_driftvector(sv, pt, missing)
    @test result == sv
end

@timed_testset "NBodyChargeCloud Units" begin

    cartesian_unit = CartesianPoint(0.01u"m", 0.0u"m", 0.05u"m")
    cartesian_float = CartesianPoint(0.01, 0.0, 0.05)
    cylindrical_unit = CylindricalPoint(0.01u"m", (π/4)u"rad", 0.05u"m")
    cylindrical_float = CylindricalPoint(0.01, π/4, 0.05)
    Edep = 1u"MeV"
    
    # Variant 1: CartesianPoint with Units, energy
    nb1 = NBodyChargeCloud(cartesian_unit, Edep)
    @test all(x -> x isa CartesianPoint, nb1.locations)
    @test isapprox(sum(nb1.energies), ustrip.(u"eV", Edep))

    # Variant 2: CylindricalPoint with Units, energy
    nb2 = NBodyChargeCloud(cylindrical_unit, Edep)
    @test all(x -> x isa CartesianPoint, nb2.locations)
    @test isapprox(sum(nb2.energies), ustrip.(u"eV", Edep))

    # Variant 3: CartesianPoint with Units, energy, N
    nb3 = NBodyChargeCloud(cartesian_unit, Edep, 10, radius = 0.001, number_of_shells = 1)
    @test all(x -> x isa CartesianPoint, nb3.locations)
    @test length(nb3.locations) >= 10
    @test isapprox(sum(nb3.energies), ustrip.(u"eV", Edep))

    # Variant 4: CylindricalPoint with Units, energy, N
    nb4 = NBodyChargeCloud(cylindrical_unit, Edep, 10, radius = 0.001, number_of_shells = 1)
    @test all(x -> x isa CartesianPoint, nb4.locations)
    @test length(nb4.locations) >= 10
    @test isapprox(sum(nb4.energies), ustrip.(u"eV", Edep))

    # Edge case, zero radius
    nb_zero = NBodyChargeCloud(cartesian_float, Edep, 5, radius = 0.0, number_of_shells = 1)
    @test all(x -> x isa CartesianPoint, nb_zero.locations)
    @test all(x -> x == nb_zero.locations[1], nb_zero.locations)  # all points coincide at the center
    @test isapprox(sum(nb_zero.energies), ustrip.(u"eV", Edep))
    @test all(x -> all(isfinite, (x.x, x.y, x.z)), nb_zero.locations)

    nb_zero_cyl = NBodyChargeCloud(cylindrical_float, Edep, 5, radius = 0.0, number_of_shells = 1)
    @test all(x -> x isa CartesianPoint, nb_zero_cyl.locations)
    @test all(x -> x == nb_zero_cyl.locations[1], nb_zero_cyl.locations) 
    @test isapprox(sum(nb_zero_cyl.energies), ustrip.(u"eV", Edep))
    @test all(x -> all(isfinite, (x.x, x.y, x.z)), nb_zero_cyl.locations)
end
