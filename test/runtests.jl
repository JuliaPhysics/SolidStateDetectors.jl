# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).
using Test

using SolidStateDetectors

using SpecialFunctions
using StaticArrays
using Statistics
using StatsBase
using Tables, TypedTables
using Unitful

include("test_utils.jl")

@timed_testset "Comparison to analytic solutions" begin
    include("test_analytic_solutions.jl")
end

@timed_testset "SOR GPU Backend" begin
    include("test_gpu_backend.jl")
end

@timed_testset "ConstructiveSolidGeometry" begin
    include("ConstructiveSolidGeometry/CSG_IO.jl")
    include("ConstructiveSolidGeometry/CSG_primitives.jl")
end

@timed_testset "Test real detectors" begin
    include("test_real_detectors.jl")
end

@timed_testset "Test charge drift" begin
    include("test_charge_drift_models.jl")
end

@timed_testset "Fano factor" begin
    

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
    include("test_io.jl")
end 

@timed_testset "Geant4 extension" begin
    if Sys.WORD_SIZE == 64 include("test_geant4.jl") end
end

display(testtimer())
