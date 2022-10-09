# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).
using Test

# using CUDAKernels # Uncomment this line in order to run all tests on (CUDA) GPU
using SolidStateDetectors

using SpecialFunctions
using StaticArrays
using StatsBase
using Tables, TypedTables
using Unitful

T = Float32
device_array_type = (@isdefined CUDAKernels) ? CUDAKernels.CUDA.CuArray : Array

@testset "Comparison to analytic solutions" begin
    include("comparison_to_analytic_solutions.jl")
end

@testset "SOR GPU Backend" begin
    include("SOR_GPU_Backend.jl")
end

@testset "ConstructiveSolidGeometry" begin
    include("ConstructiveSolidGeometry/CSG_IO.jl")
    include("ConstructiveSolidGeometry/CSG_primitives.jl")
end

T = Float32

@testset "Test real detectors" begin
    @testset "Simulate example detector: Inverted Coax" begin
        sim = Simulation{T}(SSD_examples[:InvertedCoax])
        simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 10e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
        nt = NamedTuple(sim)
        @test sim == Simulation(nt)
        deplV = estimate_depletion_voltage(sim, (verbose = false,))
        @test isapprox(deplV, 1870*u"V", atol = 5.0*u"V") 
        id = SolidStateDetectors.determine_bias_voltage_contact_id(sim.detector)
        # Check wether detector is undepleted, 10V below the previously calculated depletion voltage
        sim.detector = SolidStateDetector(sim.detector, contact_id = id, contact_potential = ustrip(deplV - deplV*0.005))
        calculate_electric_potential!(sim, depletion_handling = true)
        @test !is_depleted(sim.point_types)
        # Check wether detector is depleted, 10V above the previously calculated depletion voltage
        sim.detector = SolidStateDetector(sim.detector, contact_id = id, contact_potential = ustrip(deplV + deplV*0.005))
        calculate_electric_potential!(sim, depletion_handling = true)
        @test is_depleted(sim.point_types)
    end
    @testset "Simulate example detector: Inverted Coax (in cryostat)" begin
        sim = Simulation{T}(SSD_examples[:InvertedCoaxInCryostat])
        simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 10e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
    # @testset "Simulate example detector: Coax" begin
    #     sim = Simulation{T}(SSD_examples[:Coax])
    #     calculate_electric_potential!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
    #     calculate_electric_field!(sim)
    #     calculate_weighting_potential!(sim, 1, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
    #     calculate_weighting_potential!(sim, 2, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
    #     for i in 3:19
    #         calculate_weighting_potential!(sim, i, convergence_limit = 5e-6, refinement_limits = [0.2], verbose = false)
    #     end
    #     evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(30), 12e-3 )]))
    #     simulate!(evt, sim, Δt = 5e-10, max_nsteps = 10000)
    #     signalsum = T(0)
    #     for i in 1:length(evt.waveforms)
    #         signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
    #     end
    #     signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
    #     @info signalsum
    #     @test isapprox( signalsum, T(2), atol = 5e-3 )
    # end
    @testset "Simulate example detector: BEGe" begin
        sim = Simulation{T}(SSD_examples[:BEGe])
        calculate_electric_potential!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        calculate_electric_field!(sim)
        calculate_weighting_potential!(sim, 1, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        for i in 2:5
            calculate_weighting_potential!(sim, i, convergence_limit = 5e-6, device_array_type = device_array_type, refinement_limits = [0.2], verbose = false)
        end
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 20e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
        nt = NamedTuple(sim)
        @test sim == Simulation(nt)
        # test handling at r = 0 (see PR #322)
        calculate_weighting_potential!(sim, 4, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1, 0.05], verbose = false)
        let wp = sim.weighting_potentials[4]
            φidx1 = SolidStateDetectors.searchsortednearest(wp.grid.φ, T(deg2rad(30)))
            φidx2 = SolidStateDetectors.searchsortednearest(wp.grid.φ, T(deg2rad(210)))
            rmin = SolidStateDetectors.searchsortednearest(wp.grid.r, T(0.01))
            rmax = SolidStateDetectors.searchsortednearest(wp.grid.r, T(0.035))
            zidx = SolidStateDetectors.searchsortednearest(wp.grid.z, T(0.02))
            l = vcat(wp.data[rmin:-1:2,φidx2,zidx], wp.data[1:rmax,φidx1,zidx])
            @test all(diff(l) .< 0) # the weighting potential should monotonously decrease from rmin to rmax
        end
        @test isapprox(estimate_depletion_voltage(sim, (verbose = false,)), 0u"V", atol = 0.2u"V") # This detector has no impurity profile
    end
    @testset "Simulate example detector: HexagonalPrism" begin
        sim = Simulation{T}(SSD_examples[:Hexagon])
        simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event([CartesianPoint{T}(0, 5e-4, 1e-3)])
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
        @test isapprox(estimate_depletion_voltage(sim, (verbose = false,)), T(-13.15)*u"V", atol = 1.0u"V") 
    end
    @testset "Simulate example detector: CGD" begin
        sim = Simulation{T}(SSD_examples[:CGD])
        simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event([CartesianPoint{T}(0,2e-3,0)])
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
        nt = NamedTuple(sim)
        @test sim == Simulation(nt)
    end
    @testset "Simulate example detector: Spherical" begin
        sim = Simulation{T}(SSD_examples[:Spherical])
        simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event([CartesianPoint{T}(0,0,0)])
        simulate!(evt, sim)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end 
    @testset "Simulate example detector: Toroidal" begin
        sim = Simulation{T}(SSD_examples[:CoaxialTorus])
        SolidStateDetectors.apply_initial_state!(sim, ElectricPotential)
        simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1, 0.05, 0.02, 0.01], 
            max_tick_distance = 0.5u"mm", verbose = false)
        evt = Event([CartesianPoint{T}(0.0075,0,0)])
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
    @testset "Simulate example detector: Coaxial for partial phi range" begin
        sim = Simulation{T}(SSD_examples[:Cone2D])
        calculate_electric_potential!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = missing, verbose = false)
        
        sim_alt = Simulation{T}(SSD_examples[:Cone2D])
        calculate_electric_potential!(sim_alt, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = missing, verbose = false)
        
        # Test equal initial grid spacing in r and z, no matter the phi range
        @test sim.electric_potential.grid.r == sim_alt.electric_potential.grid.r 
        @test sim.electric_potential.grid.z == sim_alt.electric_potential.grid.z
        
        idx = findall(pt -> SolidStateDetectors.is_pn_junction_point_type(pt), sim.point_types.data)
        @test maximum(abs.(sim_alt.electric_potential.data[idx] .- sim.electric_potential.data[idx])) .< T(0.2)
    end
    @testset "Simulate example detector: SigGen PPC" begin
        sim = Simulation{T}(SSD_examples[:SigGen])
        simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 10e-3 )]))
        simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
end

@testset "Diffusion and Self-Repulsion" begin
    sim = Simulation{T}(SSD_examples[:InvertedCoax])
    simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1, 0.05], verbose = false)

    pos = CartesianPoint{T}(0.02,0,0.05); Edep = 1u"eV"
    nbcc = NBodyChargeCloud(pos, Edep, 40, radius = T(0.0005), number_of_shells = 2)

    evt = Event(nbcc)
    simulate!(evt, sim, self_repulsion = true, diffusion = true)
    signalsum = T(0)
    for i in 1:length(evt.waveforms)
        signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
    end
    signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
    @info signalsum
    @test isapprox( signalsum, T(2), atol = 5e-3 )
end

@testset "Table Simulation" begin 
    sim = Simulation{T}(SSD_examples[:InvertedCoax])
    simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)

    evt_table = Table(
        evtno = Int32[1], 
        detno = Int32[1],
        thit = [T[0] * u"s"],
        edep = [T[1] * u"eV"],
        pos = [[SVector{3, T}.(0.01, 0.01, 0.01) * u"m"]]
    )
    contact_charge_signals = simulate_waveforms(      
        evt_table,
        sim,
        max_nsteps = 4000, 
        Δt = 1u"ns", 
        number_of_carriers = 20,
        number_of_shells = 2,
        verbose = false);
    signalsum = T(0)
    for i in 1:length(contact_charge_signals.waveform)
        signalsum += abs(ustrip(contact_charge_signals.waveform[i].signal[end]))
    end
    signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
    @test isapprox( signalsum, T(2), atol = 5e-3 )
end

@testset "Isochrone" begin
    sim = Simulation{T}(SSD_examples[:IsochroneTest])
    simulate!(sim, device_array_type = device_array_type)
    calculate_electric_potential!(sim, device_array_type = device_array_type, refinement_limits = [0.2, 0.1, 0.05, 0.01])
    calculate_electric_field!(sim, n_points_in_φ = 72)
    sim.detector = SolidStateDetector(sim.detector, ADLChargeDriftModel());

    spawn_positions = CartesianPoint{T}[]
    idx_spawn_positions = CartesianIndex[]
    x_axis = T.(0:0.0005:0.029)
    z_axis = T.(0.0005:0.0005:0.029)
    for (i,x) in enumerate(x_axis)
        for (k,z) in enumerate(z_axis)
            push!(spawn_positions, CartesianPoint{T}(x,0,z))
            push!(idx_spawn_positions, CartesianIndex(i,k))
        end
    end
    length(spawn_positions)
    in_idx = findall(x -> x in sim.detector && !in(x, sim.detector.contacts), spawn_positions);
    ev = Event(spawn_positions[in_idx]);
    time_step = T(2)u"ns"
    max_nsteps = 10000
    drift_charges!(ev, sim, Δt = time_step, max_nsteps = max_nsteps, verbose = false)

    DT_h = Array{typeof(time_step),2}(fill(typeof(time_step)(NaN),length(x_axis),length(z_axis)))
    DT_e = Array{typeof(time_step),2}(fill(typeof(time_step)(NaN),length(x_axis),length(z_axis)))
    for (i, idx) in enumerate(idx_spawn_positions[in_idx])
        DT_h[idx] = length(ev.drift_paths[i].h_path)*time_step
        DT_e[idx] = length(ev.drift_paths[i].e_path)*time_step
    end

    # All charge drifts should end at some point
    @test maximum(filter(x -> !isnan(x), DT_h)) <  max_nsteps * time_step
    @test maximum(filter(x -> !isnan(x), DT_e)) <  max_nsteps * time_step
    
    # The points around the core contact should have small hole drift times
    @test maximum(filter(x -> !isnan(x), DT_h[1:20,1:20])) < 200u"ns"
    
    # The points in the upper and outer corners should have larger hole drift times
    @test minimum(filter(x -> !isnan(x), DT_h[end-14:end, end-14:end])) > 1000u"ns"
    
    # Test if all holes made it to the point contact and all electrons made it to the mantle contact
    @test all(broadcast(dp -> dp.h_path[end] in sim.detector.contacts[1], ev.drift_paths))
    @test all(broadcast(dp -> dp.e_path[end] in sim.detector.contacts[2], ev.drift_paths))
end

@testset "ADLChargeDriftModel" begin
    include("ADLChargeDriftModel.jl")
end

@testset "IO" begin
    sim = Simulation(SSD_examples[:InvertedCoax])
    
    calculate_electric_potential!(sim, verbose = false, device_array_type = device_array_type)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)
    
    calculate_electric_field!(sim)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    calculate_weighting_potential!(sim, 1, verbose = false, device_array_type = device_array_type)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    for i in findall(ismissing.(sim.weighting_potentials)) calculate_weighting_potential!(sim, i, verbose = false, device_array_type = device_array_type) end
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)
end 


