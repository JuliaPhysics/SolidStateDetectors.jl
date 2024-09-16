using OrderedCollections
using LegendHDF5IO

@testset "Read in different config file formats" begin

    # Create JSON file from YAML file
    filename = SSD_examples[:InvertedCoax]
    @test_nowarn SolidStateDetectors.yaml2json(filename)
    
    # Test read in using YAML, JSON and SigGen using Dict, OrderedDict and OrderedDict{String, Any}
    for f in (filename, replace(filename, ".yaml" => ".json"), SSD_examples[:SigGen])
        @test isfile(f)
        @testset "$(split(f, ".")[end])" begin
            @test_nowarn SolidStateDetectors.parse_config_file(f)
            @test_nowarn SolidStateDetectors.parse_config_file(f, dicttype = OrderedDict)
            @test_nowarn SolidStateDetectors.parse_config_file(f, dicttype = OrderedDict{String, Any})
        end
    end
end

sim = Simulation(SSD_examples[:InvertedCoax])

@testset "Conversion Simulation <--> NamedTuple" begin
    timed_calculate_electric_potential!(sim, verbose = false)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    timed_calculate_electric_field!(sim)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    timed_calculate_weighting_potential!(sim, 1, verbose = false)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    for i in findall(ismissing.(sim.weighting_potentials)) timed_calculate_weighting_potential!(sim, i, verbose = false) end
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)
end

@testset "Test LegendHDF5IO" begin
    @test_nowarn ssd_write("legendhdf5io_test.lh5", sim)
    @test_logs (:warn, "Destination `legendhdf5io_test.lh5` already exists. Overwriting...") ssd_write("legendhdf5io_test.lh5", sim)
    @test sim == ssd_read("legendhdf5io_test.lh5", Simulation)
end