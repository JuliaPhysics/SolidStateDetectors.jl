# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).
using Test

# using CUDAKernels # Uncomment this line in order to run all tests on (CUDA) GPU
using SolidStateDetectors

using SpecialFunctions
using StaticArrays
using Statistics
using StatsBase
using Tables, TypedTables
using Unitful
using TimerOutputs


T = Float32
device_array_type = (@isdefined CUDAKernels) ? CUDAKernels.CUDA.CuArray : Array


testtimer() = get_timer("_default_testtimer_")

macro timed_testset(title, body)
    quote
        tmr = testtimer()
        _title = $(esc(title))
        @timeit tmr "$_title" begin
            @testset "$_title" begin
                $(esc(body))
            end
        end
    end
end


function timed_calculate_electric_potential!(args...; kwargs...)
    @timeit testtimer() "calculate_electric_potential!" begin
        calculate_electric_potential!(args...; kwargs...)
    end
end

function timed_calculate_electric_field!(args...; kwargs...)
    @timeit testtimer() "calculate_electric_field!" begin
        calculate_electric_field!(args...; kwargs...)
    end
end

function timed_calculate_weighting_potential!(args...; kwargs...)
    @timeit testtimer() "calculate_weighting_potential!" begin
        calculate_weighting_potential!(args...; kwargs...)
    end
end

function timed_simulate_waveforms(args...; kwargs...)
    @timeit testtimer() "simulate_waveforms" begin
        simulate_waveforms(args...; kwargs...)
    end
end

function timed_simulate!(args...; kwargs...)
    @timeit testtimer() "simulate!" begin
        simulate!(args...; kwargs...)
    end
end

function timed_estimate_depletion_voltage(args...; kwargs...)
    @timeit testtimer() "estimate_depletion_voltage" begin
        estimate_depletion_voltage(args...; kwargs...)
    end
end


@timed_testset "Comparison to analytic solutions" begin
    include("comparison_to_analytic_solutions.jl")
end

@timed_testset "SOR GPU Backend" begin
    include("SOR_GPU_Backend.jl")
end

@timed_testset "ConstructiveSolidGeometry" begin
    include("ConstructiveSolidGeometry/CSG_IO.jl")
    include("ConstructiveSolidGeometry/CSG_primitives.jl")
end

T = Float32

@timed_testset "Test real detectors" begin
    @timed_testset "Simulate example detector: Inverted Coax" begin
        sim = Simulation{T}(SSD_examples[:InvertedCoax])
        timed_simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1, 0.05, 0.03, 0.02, 0.01], verbose = false)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 10e-3 )]))
        timed_simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
        nt = NamedTuple(sim)
        @test sim == Simulation(nt)
        deplV = timed_estimate_depletion_voltage(sim, verbose = false)
        @test isapprox(deplV, 1870*u"V", atol = 5.0*u"V") 
        id = SolidStateDetectors.determine_bias_voltage_contact_id(sim.detector)
        # Check wether detector is undepleted, 10V below the previously calculated depletion voltage
        sim.detector = SolidStateDetector(sim.detector, contact_id = id, contact_potential = ustrip(deplV - deplV*0.005))
        timed_calculate_electric_potential!(sim, depletion_handling = true)
        @test !is_depleted(sim.point_types)
        # Check wether detector is depleted, 10V above the previously calculated depletion voltage
        sim.detector = SolidStateDetector(sim.detector, contact_id = id, contact_potential = ustrip(deplV + deplV*0.005))
        timed_calculate_electric_potential!(sim, depletion_handling = true)
        @test is_depleted(sim.point_types)
    end
    @timed_testset "Simulate example detector: Inverted Coax (in cryostat)" begin
        sim = Simulation{T}(SSD_examples[:InvertedCoaxInCryostat])
        timed_simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 10e-3 )]))
        timed_simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
    # @timed_testset "Simulate example detector: Coax" begin
    #     sim = Simulation{T}(SSD_examples[:Coax])
    #     timed_calculate_electric_potential!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
    #     timed_calculate_electric_field!(sim)
    #     timed_calculate_weighting_potential!(sim, 1, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
    #     timed_calculate_weighting_potential!(sim, 2, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1], verbose = false)
    #     for i in 3:19
    #         timed_calculate_weighting_potential!(sim, i, convergence_limit = 5e-6, refinement_limits = [0.2], verbose = false)
    #     end
    #     evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(30), 12e-3 )]))
    #     timed_simulate!(evt, sim, Δt = 5e-10, max_nsteps = 10000)
    #     signalsum = T(0)
    #     for i in 1:length(evt.waveforms)
    #         signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
    #     end
    #     signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
    #     @info signalsum
    #     @test isapprox( signalsum, T(2), atol = 5e-3 )
    # end
    @timed_testset "Simulate example detector: BEGe" begin
        sim = Simulation{T}(SSD_examples[:BEGe])
        timed_calculate_electric_potential!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        timed_calculate_electric_field!(sim)
        timed_calculate_weighting_potential!(sim, 1, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        for i in 2:5
            timed_calculate_weighting_potential!(sim, i, convergence_limit = 5e-6, device_array_type = device_array_type, refinement_limits = [0.2], verbose = false)
        end
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 20e-3 )]))
        timed_simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
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
        timed_calculate_weighting_potential!(sim, 4, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1, 0.05], verbose = false)
        @timed_testset "monotonous decrease of WP" let wp = sim.weighting_potentials[4]
            φidx1 = SolidStateDetectors.searchsortednearest(wp.grid.φ, T(deg2rad(30)))
            φidx2 = SolidStateDetectors.searchsortednearest(wp.grid.φ, T(deg2rad(210)))
            rmin = SolidStateDetectors.searchsortednearest(wp.grid.r, T(0.01))
            rmax = SolidStateDetectors.searchsortednearest(wp.grid.r, T(0.035))
            zidx = SolidStateDetectors.searchsortednearest(wp.grid.z, T(0.02))
            l = vcat(wp.data[rmin:-1:2,φidx2,zidx], wp.data[1:rmax,φidx1,zidx])
            @test all(diff(l) .< 0) # the weighting potential should monotonously decrease from rmin to rmax
        end
        @test isapprox(timed_estimate_depletion_voltage(sim, verbose = false), 0u"V", atol = 0.2u"V") # This detector has no impurity profile
    end
    @timed_testset "Simulate example detector: HexagonalPrism" begin
        sim = Simulation{T}(SSD_examples[:Hexagon])
        timed_simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1, 0.05, 0.02], verbose = false)
        evt = Event([CartesianPoint{T}(0, 5e-4, 1e-3)])
        timed_simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
        @test isapprox(timed_estimate_depletion_voltage(sim, verbose = false), T(-13.15)*u"V", atol = 1.0u"V") 
    end
    @timed_testset "Simulate example detector: CGD" begin
        sim = Simulation{T}(SSD_examples[:CGD])
        timed_simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event([CartesianPoint{T}(0,2e-3,0)])
        timed_simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
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
    @timed_testset "Simulate example detector: Spherical" begin
        sim = Simulation{T}(SSD_examples[:Spherical])
        timed_simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event([CartesianPoint{T}(0,0,0)])
        timed_simulate!(evt, sim)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end 
    @timed_testset "Simulate example detector: Toroidal" begin
        sim = Simulation{T}(SSD_examples[:CoaxialTorus])
        SolidStateDetectors.apply_initial_state!(sim, ElectricPotential)
        timed_simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1, 0.05, 0.02, 0.01], 
            max_tick_distance = 0.5u"mm", verbose = false)
        evt = Event([CartesianPoint{T}(0.0075,0,0)])
        timed_simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
    @timed_testset "Simulate example detector: Coaxial for partial phi range" begin
        sim = Simulation{T}(SSD_examples[:Cone2D])
        timed_calculate_electric_potential!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = missing, verbose = false)
        
        sim_alt = Simulation{T}(SSD_examples[:Cone2D])
        timed_calculate_electric_potential!(sim_alt, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = missing, verbose = false)
        
        # Test equal initial grid spacing in r and z, no matter the phi range
        @test sim.electric_potential.grid.r == sim_alt.electric_potential.grid.r 
        @test sim.electric_potential.grid.z == sim_alt.electric_potential.grid.z
        
        idx = findall(pt -> SolidStateDetectors.is_pn_junction_point_type(pt), sim.point_types.data)
        @test maximum(abs.(sim_alt.electric_potential.data[idx] .- sim.electric_potential.data[idx])) .< T(0.2)
    end
    @timed_testset "Simulate example detector: SigGen PPC" begin
        sim = Simulation{T}(SSD_examples[:SigGen])
        timed_simulate!(sim, convergence_limit = 1e-6, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)
        evt = Event(CartesianPoint.([CylindricalPoint{T}(20e-3, deg2rad(10), 10e-3 )]))
        timed_simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
        signalsum = T(0)
        for i in 1:length(evt.waveforms)
            signalsum += abs(ustrip(evt.waveforms[i].signal[end]))
        end
        signalsum *= inv(ustrip(SolidStateDetectors._convert_internal_energy_to_external_charge(sim.detector.semiconductor.material)))
        @info signalsum
        @test isapprox( signalsum, T(2), atol = 5e-3 )
    end
end

@timed_testset "Charge Drift" begin
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
        @test signalsum < T(2)
    end

    @timed_testset "ADLChargeDriftModel" begin
        include("ADLChargeDriftModel.jl")
    end
end

@timed_testset "Fano factor" begin
    

end

@timed_testset "Table Simulation" begin 
    sim = Simulation{T}(SSD_examples[:InvertedCoax])
    timed_simulate!(sim, convergence_limit = 1e-5, device_array_type = device_array_type, refinement_limits = [0.2, 0.1], verbose = false)

    evt_table = Table(
        evtno = Int32[1], 
        detno = Int32[1],
        thit = [T[0] * u"s"],
        edep = [T[1] * u"eV"],
        pos = [[SVector{3, T}.(0.01, 0.01, 0.01) * u"m"]]
    )
    contact_charge_signals = timed_simulate_waveforms(      
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

@timed_testset "Isochrone" begin
    T = Float64 # to prevent rounding errors that might cause the final test to fail
    sim = Simulation{T}(SSD_examples[:IsochroneTest])
    timed_simulate!(sim, device_array_type = device_array_type)
    timed_calculate_electric_potential!(sim, device_array_type = device_array_type, refinement_limits = [0.2, 0.1, 0.05, 0.01])
    timed_calculate_electric_field!(sim, n_points_in_φ = 72)
    sim.detector = SolidStateDetector(sim.detector, ADLChargeDriftModel{T}());

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

@timed_testset "IO" begin
    include("IO.jl")
end 

@timed_testset "Geant4 extension" begin
    if Sys.WORD_SIZE == 64 include("Geant4.jl") end
end

display(testtimer())
