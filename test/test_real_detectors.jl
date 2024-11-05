# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors
using Unitful

# include("test_utils.jl")

T = Float32

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

