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

@testset "Surface refinement edge cases" begin
    T = Float64
    # Surface model but no spacing_surface_refinement
    begin
        sim = @test_nowarn Simulation{T}(SSD_examples[:TrueCoaxial])
        
        @test_logs (:warn, r"Surface model detected but `spacing_surface_refinement` is not defined") begin
            timed_calculate_electric_potential!(sim, verbose = false, depletion_handling = true)
        end
    end
    
    # Wrong refinement limits
    begin
        sim_cyl = @test_nowarn Simulation{T}(SSD_examples[:IVCIlayer])
        sim_cyl.spacing_surface_refinement = (1e-3, 1e-3, 1e-3)
        # Test normal behaviour Cylindrical
        timed_calculate_electric_potential!(sim_cyl, verbose = false, depletion_handling = true)
        grid_ax1_cyl = length(sim_cyl.electric_potential.grid[1])
        grid_ax2_cyl = length(sim_cyl.electric_potential.grid[2])
        grid_ax3_cyl = length(sim_cyl.electric_potential.grid[3])
        
        # --- Case: < 3 refinements ---
        bad_limits = [0.2, 0.1]
        @test_logs (:warn, r"Surface model detected") (:warn, r"Falling back to default") match_mode=:any begin
            timed_calculate_electric_potential!(sim_cyl, refinement_limits=bad_limits, verbose = false, depletion_handling = true)
        end
        
        # --- Case: last refinement > 0.05 ---
        bad_limits = [0.2, 0.1, 0.09]
        @test_logs (:warn, r"Surface model detected:") (:warn, r"Falling back to default") match_mode=:any begin
            timed_calculate_electric_potential!(sim_cyl, refinement_limits=bad_limits, verbose = false, depletion_handling = true)
        end

        # Test normal behaviour Cartesian
        config_dict = SolidStateDetectors.parse_config_file(SSD_examples[:IVCIlayer])
        config_dict["grid"]["coordinates"] = "cartesian"
        config_dict["grid"]["axes"]["x"] = Dict(
            "from" => "-40",
            "to"   => "40",
            "boundaries"   => "inf"
        )
        config_dict["grid"]["axes"]["y"] = Dict(
            "from" => "-40",
            "to"   => "40",
            "boundaries"   => "inf"
        )
        config_dict["grid"]["axes"]["z"] = Dict(
            "from" => "-10",
            "to"   => "90",
            "boundaries"   => "inf"
        )
        sim_cart = @test_nowarn Simulation{T}(config_dict)
        sim_cart.spacing_surface_refinement = (1e-3, 1e-3, 1e-3)
        timed_calculate_electric_potential!(sim_cart, verbose = false, depletion_handling = true)
        grid_ax1_cart = length(sim_cart.electric_potential.grid[1])
        grid_ax2_cart = length(sim_cart.electric_potential.grid[2])
	grid_ax3_cart = length(sim_cart.electric_potential.grid[3])
        
        @test grid_ax1_cyl == 74
        @test grid_ax2_cyl == 1
        @test grid_ax3_cyl == 143
        @test grid_ax1_cart == 178
        @test grid_ax2_cart == 178
        @test grid_ax3_cart == 149
    end
    
    # Spacing out of bounds
    sim = @test_nowarn Simulation{T}(SSD_examples[:IVCIlayer])
    large_spacing = (1e-2, 1e-2, 1e-2)
    sim.spacing_surface_refinement = large_spacing
    
    @test_logs (
        (:warn, r"out of bounds"s)
    ) match_mode=:any begin
        timed_calculate_electric_potential!(sim, verbose=false, depletion_handling=true)
    end
          
    # Cylindrical grid: φ-axis not refined
    begin
        config_dict = SolidStateDetectors.parse_config_file(SSD_examples[:IVCIlayer])
        
        config_dict["grid"]["axes"]["phi"] = Dict(
            "from" => "0",
            "to"   => "120"
        )
        
        sim_phi = @test_nowarn Simulation{T}(config_dict)
        sim_phi.spacing_surface_refinement = nothing
        timed_calculate_electric_potential!(sim_phi, verbose = false, depletion_handling = true)
        notref_phi_len = length(sim_phi.electric_potential.grid.axes[2].ticks)

        sim_phi.spacing_surface_refinement = (1e-3, 1e-3, 1e-3)
        timed_calculate_electric_potential!(sim_phi, verbose = false, depletion_handling = true)
        ref_phi_len = length(sim_phi.electric_potential.grid.axes[2].ticks)
        
        @test ref_phi_len == notref_phi_len
    end
end
