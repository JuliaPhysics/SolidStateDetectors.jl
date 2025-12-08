using SolidStateDetectors
using Test

@timed_testset "Test simulating in two refinement steps" begin

    T = Float32
    sim = Simulation{T}("test_config_files/BEGe_01.yaml")
    timed_calculate_electric_potential!(sim, refinement_limits = [0.2, 0.1], depletion_handling = true)

    sim2 = Simulation{T}("test_config_files/BEGe_01.yaml")
    timed_calculate_electric_potential!(sim2, refinement_limits = [0.2], depletion_handling = true)
    timed_calculate_electric_potential!(sim2, refinement_limits = [0.1], depletion_handling = true, initialize = false)

    @test size(sim.electric_potential.data) == size(sim2.electric_potential.data)
    @test all(isapprox.(sim.electric_potential.data, sim2.electric_potential.data))
    @test sim.electric_potential.grid == sim2.electric_potential.grid

    timed_calculate_weighting_potential!(sim, 1, refinement_limits = [0.2, 0.1], depletion_handling = true)
    timed_calculate_weighting_potential!(sim2, 1, refinement_limits = [0.2], depletion_handling = true)
    timed_calculate_weighting_potential!(sim2, 1, refinement_limits = [0.1], depletion_handling = true, initialize = false)

    @test size(sim.weighting_potentials[1].data) == size(sim2.weighting_potentials[1].data)
end

@timed_testset "Test compute min_tick_distance" begin

    CS = SolidStateDetectors.Cylindrical
    T = Float64
    
    #Test for compute_min_tick_distance
    
    struct TestGrid
           axes::NTuple{3, Tuple{T, T}}
    end
    
    grid_small = TestGrid((
        (0.0, 0.01),       # r: 0 → 10 mm
        (0.0, 120),        # phi
        (-0.0005, 0.0005)  # z: -500 µm → 500 µm
              ))
    min_tick_small = SolidStateDetectors.compute_min_tick_distance(grid_small, T, CS)
    @test min_tick_small[3] == 1.0e-6
    
    grid_large = TestGrid((
        (0.0, 0.1),        # r: 0 → 100 mm
        (0.0, 120),        # phi
        (0.0, 0.1)         # z: 0 → 100 mm
    ))
    min_tick_large = SolidStateDetectors.compute_min_tick_distance(grid_large, T, CS)
    @test min_tick_large[3] == 1.0e-5
    
    #Test for compute_min_tick_distance in simulation 

    sim = Simulation{T}("test_config_files/BEGe_03.yaml")
    timed_calculate_electric_potential!(sim)

    @test length(sim.electric_potential.grid[1]) == 38
    @test length(sim.electric_potential.grid[2]) == 6
    @test length(sim.electric_potential.grid[3]) == 68
end
