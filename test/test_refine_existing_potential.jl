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

    T = Float64
    
    #Test for compute_min_tick_distance
    axr = SolidStateDetectors.DiscreteAxis(
        0.0, 0.01, :left, :right, :L, :R, collect(range(0.0, stop=0.01, length=10)))              # r: 0 → 10 mm
    axφ = SolidStateDetectors.DiscreteAxis(
        0.0, 2.094, :left, :right, :L, :R, collect(range(0.0, stop=2.094, length=6)))             # φ: 0 → 120° in radians
    axz = SolidStateDetectors.DiscreteAxis(
        -0.0005, 0.0005, :left, :right, :L, :R, collect(range(-0.0005, stop=0.0005, length=20)))  # z: -500 µm → 500 µm
    
    grid_small = SolidStateDetectors.CylindricalGrid{T}((axr, axφ, axz))
    min_tick_small = SolidStateDetectors.compute_min_tick_distance(grid_small)
    @test min_tick_small[3] ≈ 1e-6

    axr_large = SolidStateDetectors.DiscreteAxis(
        0.0, 0.1, :left, :right, :L, :R, collect(range(0.0, stop=0.1, length=10)))                # r: 0 → 100 mm
    axφ_large = SolidStateDetectors.DiscreteAxis(  
        0.0, 2.094, :left, :right, :L, :R, collect( range(0.0, stop=2.094, length=6)))            # φ: 0 → 120° in radians
    axz_large = SolidStateDetectors.DiscreteAxis(  
        0.0, 0.1, :left, :right, :L, :R, collect(range(0.0, stop=0.1, length=10)))                # z: 0 → 100 mm
    
    grid_large = SolidStateDetectors.CylindricalGrid{T}((axr_large, axφ_large, axz_large))
    min_tick_large = SolidStateDetectors.compute_min_tick_distance(grid_large)
    @test min_tick_large[3] == 1.0e-5

    axx = SolidStateDetectors.DiscreteAxis(
        0.0, 0.01, :left, :right, :L, :R, collect(range(0.0, stop=0.01, length=10)) )              # x: 0 → 10 mm
    axy = SolidStateDetectors.DiscreteAxis(
        0.0, 0.01, :left, :right, :L, :R, collect(range(0.0, stop=0.01, length=10)) )              # y: 0 → 10 mm
    axz = SolidStateDetectors.DiscreteAxis(
        -0.0005, 0.0005, :left, :right, :L, :R, collect(range(-0.0005, stop=0.0005, length=20)) )  # z: -500 µm → 500 µm
    
    grid_small = SolidStateDetectors.CartesianGrid3D{T}((axx, axy, axz))
    min_tick_small = SolidStateDetectors.compute_min_tick_distance(grid_small)
    @test min_tick_small[3] ≈ 1e-6
    
    axx_large = SolidStateDetectors.DiscreteAxis(
        0.0, 0.1, :left, :right, :L, :R, collect(range(0.0, stop=0.1, length=10)) )                # x: 0 → 100 mm
    axy_large = SolidStateDetectors.DiscreteAxis(
        0.0, 0.1, :left, :right, :L, :R, collect(range(0.0, stop=0.1, length=10)) )                # y: 0 → 100 mm
    axz_large = SolidStateDetectors.DiscreteAxis(
        0.0, 0.1, :left, :right, :L, :R, collect(range(0.0, stop=0.1, length=10)) )                # z: 0 → 100 mm
    
    grid_large = SolidStateDetectors.CartesianGrid3D{T}((axx_large, axy_large, axz_large))
    min_tick_large = SolidStateDetectors.compute_min_tick_distance(grid_large)
    @test min_tick_large[3] == 1.0e-5
    
    #Test for compute_min_tick_distance in simulation 
    # Cylindrical Case
    cf = SolidStateDetectors.parse_config_file("test_config_files/BEGe_01.yaml")
    
    cf["units"]["length"] = "µm"
    cf["detectors"][1]["semiconductor"]["impurity_density"] = Dict(
        "name" => "constant",
        "value" => 0
    )
    
    sim = Simulation{T}(cf)
    timed_calculate_electric_potential!(sim)
    @test length(sim.electric_potential.grid[1]) == 36
    @test length(sim.electric_potential.grid[2]) == 1
    @test length(sim.electric_potential.grid[3]) == 62
    
    # Cartesian Case um
    cf_cart = deepcopy(cf)
    cf_cart["grid"]["coordinates"] = "cartesian"
    cf_cart["grid"]["axes"] = Dict(
        "x" => Dict("to" => 60),
        "y" => Dict("to" => 60),
        "z" => Dict("to" => 60)
    )
    cf_cart["units"]["length"] = "µm"
    cf_cart["detectors"][1]["semiconductor"]["impurity_density"] = Dict(
        "name" => "constant",
        "value" => 0
    )
    sim_cart = Simulation{T}(cf_cart)
    timed_calculate_electric_potential!(sim_cart)
    
    @test length(sim_cart.electric_potential.grid[1]) == 60
    @test length(sim_cart.electric_potential.grid[2]) == 60
    @test length(sim_cart.electric_potential.grid[3]) == 62

    grid_lengths = length.(sim_cart.electric_potential.grid)
    @test maximum(grid_lengths) - minimum(grid_lengths) <= 2

    # Cartesian Case mm
    cf_cart_mm = deepcopy(cf_cart)
    cf_cart_mm["units"]["length"] = "mm"
    sim_cart_mm = Simulation{T}(cf_cart_mm)
    timed_calculate_electric_potential!(sim_cart_mm)

    @test length(sim_cart_mm.electric_potential.grid[1]) == 60
    @test length(sim_cart_mm.electric_potential.grid[2]) == 60
    @test length(sim_cart_mm.electric_potential.grid[3]) == 64
end
